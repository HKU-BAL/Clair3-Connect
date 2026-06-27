use crate::fa_pipeline::{connect_predict_client, predict_and_emit};
use anyhow::{Context, Result};
use clair3_connect_core::full_alignment::driver::{calculate_full_alignment_streaming, FaInputs};
use clair3_connect_core::full_alignment::matrix::FaMatrix;
use clair3_connect_core::hts::region::Region;
use clair3_connect_core::pileup::kernel::{calculate_pileup, PileupParams};
use clair3_connect_core::vcf::header::{build_vcf_header, VcfHeaderInputs};
use clair3_connect_core::vcf::output_with::OutputWithConfig;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

#[derive(Clone)]
pub struct FaOnlyInputs<'a> {

    pub bam_path: &'a Path,

    pub fasta_path: &'a Path,

    pub rois: Vec<Region>,

    pub server_endpoint: String,

    pub fa_model_id: String,

    pub sample_name: String,

    pub clair3_version: String,

    pub matrix_depth: usize,

    pub max_indel_length: usize,

    pub min_mq: u8,

    pub pileup_params: PileupParams,

    pub output_config: OutputWithConfig,

    pub batch_size: usize,

    pub fa_window_bp: u64,

    pub fa_threads: usize,

    pub aead_key: Option<clair3_connect_transport::aead::AeadKey>,

    pub tls: Option<clair3_connect_transport::tls::ClientTlsOptions>,
}

pub async fn run_fa_only_pipeline<W: Write>(
    inputs: &FaOnlyInputs<'_>,
    out: &mut W,
) -> Result<()> {

    let header = build_vcf_header(&VcfHeaderInputs {
        reference_file_path: Some(inputs.fasta_path),
        cmdline_file: None,
        sample_name: &inputs.sample_name,
        clair3_version: &inputs.clair3_version,
        gvcf: false,
    })?;
    out.write_all(header.text.as_bytes())?;
    out.write_all(b"\n")?;

    if inputs.rois.is_empty() {
        tracing::info!("fa-only: no ROIs given, header-only output");
        return Ok(());
    }

    let client = connect_predict_client(
        &inputs.server_endpoint,
        inputs.tls.as_ref(),
        inputs.aead_key.as_ref(),
    )
    .await?;

    let window_bp = if inputs.fa_window_bp == 0 {
        500_000
    } else {
        inputs.fa_window_bp
    };

    let mut contig_rank: HashMap<String, usize> = HashMap::new();
    for roi in &inputs.rois {
        let n = contig_rank.len();
        contig_rank.entry(roi.contig.clone()).or_insert(n);
    }

    let mut windows: Vec<(usize, Region)> = Vec::new();
    for roi in &inputs.rois {
        let rank = contig_rank[&roi.contig];
        let mut ws = roi.start;
        while ws < roi.end {
            let we = (ws + window_bp).min(roi.end);
            windows.push((
                rank,
                Region {
                    contig: roi.contig.clone(),
                    start: ws,
                    end: we,
                },
            ));
            ws = we;
        }
    }
    let total_windows = windows.len();

    let cpu_concurrency = inputs.fa_threads.max(1);
    let cpu_sem = std::sync::Arc::new(tokio::sync::Semaphore::new(cpu_concurrency));
    let predict_sem = std::sync::Arc::new(tokio::sync::Semaphore::new(
        crate::fa_pipeline::FA_PREDICT_CONCURRENCY,
    ));

    let mut handles = Vec::with_capacity(total_windows);
    for (rank, win) in windows {

        let cpu_permit = cpu_sem
            .clone()
            .acquire_owned()
            .await
            .expect("cpu semaphore closed");
        let client = client.clone();
        let predict_sem = predict_sem.clone();
        let bam = inputs.bam_path.to_path_buf();
        let fasta = inputs.fasta_path.to_path_buf();
        let pileup_params = inputs.pileup_params.clone();
        let model_id = inputs.fa_model_id.clone();
        let output_config = inputs.output_config.clone();
        let matrix_depth = inputs.matrix_depth;
        let batch_size = inputs.batch_size;
        let min_mq = inputs.min_mq;
        let max_indel_length = inputs.max_indel_length;
        handles.push(tokio::spawn(async move {
            let _cpu_permit = cpu_permit;
            let win_build = win.clone();

            let built = tokio::task::spawn_blocking(
                move || -> Result<Option<(Vec<u64>, FaMatrix)>> {
                    let scan = calculate_pileup(&win_build, &bam, &fasta, &pileup_params)
                        .with_context(|| {
                            format!(
                                "fa-only candidate scan {}:{}-{}",
                                win_build.contig, win_build.start, win_build.end
                            )
                        })?;
                    let mut cands: Vec<u64> = scan
                        .candidates
                        .iter()
                        .map(|c| (c.pos_1based - 1) as u64)
                        .collect();
                    drop(scan);
                    cands.sort_unstable();
                    cands.dedup();
                    if cands.is_empty() {
                        return Ok(None);
                    }
                    let fa_inputs = FaInputs {
                        bam_path: &bam,
                        fasta_path: &fasta,
                        region: (
                            win_build.contig.clone(),
                            win_build.start as i64,
                            win_build.end as i64,
                        ),
                        candidates: &cands,
                        min_mq,
                        matrix_depth,
                        max_indel_length,
                    };
                    let fa = calculate_full_alignment_streaming(&fa_inputs)
                        .context("calculate_full_alignment_streaming")?;
                    Ok(Some((cands, fa)))
                },
            )
            .await
            .map_err(|e| anyhow::anyhow!("fa build join: {e}"))??;

            let (rows, ncands) = match built {
                None => (Vec::new(), 0usize),
                Some((cands, fa)) => {
                    let rows = predict_and_emit(
                        &client,
                        &fa,
                        &cands,
                        &win.contig,
                        &model_id,
                        matrix_depth,
                        batch_size,
                        &output_config,
                        &predict_sem,
                    )
                    .await?;
                    (rows, cands.len())
                }
            };
            Ok::<(usize, Vec<(u64, String)>, usize), anyhow::Error>((rank, rows, ncands))
        }));
    }

    let mut emitted: Vec<(usize, u64, String)> = Vec::new();
    let mut total_cands = 0usize;
    for h in handles {
        let (rank, rows, ncands) = h
            .await
            .map_err(|e| anyhow::anyhow!("fa-only window join: {e}"))??;
        total_cands += ncands;
        for (pos, row) in rows {
            emitted.push((rank, pos, row));
        }
    }

    emitted.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    for (_, _, row) in &emitted {
        out.write_all(row.as_bytes())?;
    }
    tracing::info!(
        "fa-only: {} ROIs, {} windows, {} candidates, {} rows emitted ({} build threads)",
        inputs.rois.len(),
        total_windows,
        total_cands,
        emitted.len(),
        cpu_concurrency
    );
    Ok(())
}
