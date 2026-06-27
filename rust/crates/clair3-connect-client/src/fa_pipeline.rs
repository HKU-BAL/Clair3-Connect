use anyhow::{anyhow, Context, Result};
use clair3_connect_core::full_alignment::driver::{
    calculate_full_alignment_no_haplotag, calculate_full_alignment_streaming,
    calculate_full_alignment_with_haplotag, FaInputs,
};
use clair3_connect_core::full_alignment::phase_block_index::Variant;
use clair3_connect_core::hts::region::Region;
use clair3_connect_core::vcf::header::{build_vcf_header, VcfHeaderInputs};
use clair3_connect_core::vcf::output_with::{output_with, OutputWithConfig, OutputWithInputs};
use clair3_connect_transport::codec::decode_f32;
use clair3_connect_transport::predict_client::{make_request_i8, PredictClient};
use std::io::Write;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct FaPipelineInputs<'a> {

    pub bam_path: &'a Path,

    pub fasta_path: &'a Path,

    pub region: Region,

    pub candidates: Vec<u64>,

    pub phased_variants: Option<Vec<Variant>>,

    pub server_endpoint: String,

    pub model_id: String,

    pub sample_name: String,

    pub clair3_version: String,

    pub matrix_depth: usize,

    pub max_indel_length: usize,

    pub min_mq: u8,

    pub output_config: OutputWithConfig,

    pub batch_size: usize,

    pub aead_key: Option<clair3_connect_transport::aead::AeadKey>,

    pub tls: Option<clair3_connect_transport::tls::ClientTlsOptions>,

    pub low_mem: bool,

    pub fa_window_bp: u64,

    pub streaming: bool,
}

pub(crate) async fn fa_window_rows(
    client: &mut PredictClient,
    candidates: &[u64],
    region: &Region,
    inputs: &FaPipelineInputs<'_>,
) -> Result<Vec<(u64, String)>> {
    let fa_inputs = FaInputs {
        bam_path: inputs.bam_path,
        fasta_path: inputs.fasta_path,
        region: (region.contig.clone(), region.start as i64, region.end as i64),
        candidates,
        min_mq: inputs.min_mq,
        matrix_depth: inputs.matrix_depth,
        max_indel_length: inputs.max_indel_length,
    };
    let fa = if inputs.streaming {

        calculate_full_alignment_streaming(&fa_inputs)
            .context("calculate_full_alignment_streaming")?
    } else if let Some(variants) = &inputs.phased_variants {
        calculate_full_alignment_with_haplotag(&fa_inputs, variants)
            .context("calculate_full_alignment_with_haplotag")?
    } else {
        calculate_full_alignment_no_haplotag(&fa_inputs)
            .context("calculate_full_alignment_no_haplotag")?
    };

    let elements_per_row: usize =
        inputs.matrix_depth * fa.no_of_positions * fa.channel_size;
    let n = fa.candidate_num;
    const OUT_COLS: usize = 90;
    const FA_PREDICT_CONCURRENCY: usize = 8;
    let mut y_all: Vec<f32> = vec![0.0; n * OUT_COLS];
    let batch_rows: Vec<usize> = (0..n).step_by(inputs.batch_size).collect();
    for wave in batch_rows.chunks(FA_PREDICT_CONCURRENCY) {
        let mut handles = Vec::with_capacity(wave.len());
        for &row in wave {
            let end = (row + inputs.batch_size).min(n);
            let batch = end - row;
            let slice = &fa.matrix[row * elements_per_row..end * elements_per_row];
            let req = make_request_i8(
                &format!("fa_batch{}", row),
                &inputs.model_id,
                vec![
                    batch as i64,
                    inputs.matrix_depth as i64,
                    fa.no_of_positions as i64,
                    fa.channel_size as i64,
                ],
                slice,
            );
            let mut c = client.clone();
            handles.push(tokio::spawn(async move {
                let resp = c
                    .predict(req)
                    .await
                    .map_err(|e| anyhow!("fa predict batch {}: {}", row, e))?;
                Ok::<(usize, usize, Vec<f32>), anyhow::Error>((
                    row,
                    batch,
                    decode_f32(&resp.y_payload),
                ))
            }));
        }
        for h in handles {
            let (row, batch, y_batch) =
                h.await.map_err(|e| anyhow!("fa predict join: {e}"))??;
            let expected = batch * OUT_COLS;
            if y_batch.len() != expected {
                return Err(anyhow!(
                    "fa Y batch {} unexpected size: got {}, want {}",
                    row,
                    y_batch.len(),
                    expected
                ));
            }
            y_all[row * OUT_COLS..row * OUT_COLS + y_batch.len()].copy_from_slice(&y_batch);
        }
    }

    let use_fa_path = inputs.output_config.add_indel_length;
    let mut emitted: Vec<(u64, String)> = Vec::with_capacity(fa.candidate_num);
    for i in 0..fa.candidate_num {
        let r = i * 90;
        let gt21 = &y_all[r..r + 21];
        let geno = &y_all[r + 21..r + 24];
        let vlen1 = &y_all[r + 24..r + 57];
        let vlen2 = &y_all[r + 57..r + 90];
        let candidate_0based = candidates[i];
        let candidate_1based = candidate_0based + 1;
        let alt_info = &fa.alt_info_strings[i];
        let mut parts = alt_info.splitn(4, '-');
        let _pos_field = parts.next().unwrap_or("0");
        let depth_field = parts.next().unwrap_or("0");
        let ref_byte = parts.next().and_then(|s| s.bytes().next()).unwrap_or(b'N');
        let payload_after = parts.next().unwrap_or("");
        let ref_char = ref_byte as char;
        let chr_pos_seq =
            format!("{}:{}:{}", inputs.region.contig, candidate_1based, ref_char);
        let alt_info_for_output = format!("{}-{}", depth_field, payload_after);
        let inputs_row = OutputWithInputs {
            chr_pos_seq: &chr_pos_seq,
            alt_info: &alt_info_for_output,
            gt21_probabilities: gt21,
            genotype_probabilities: geno,
            vlen_heads: if use_fa_path { Some((vlen1, vlen2)) } else { None },
        };
        match output_with(inputs_row, &inputs.output_config) {
            Ok(Some(row_text)) => emitted.push((candidate_1based, row_text)),
            Ok(None) => {}
            Err(e) => tracing::warn!("fa output_with row {}: {e}", i),
        }
    }
    Ok(emitted)
}

pub(crate) const FA_PREDICT_CONCURRENCY: usize = 8;

#[allow(clippy::too_many_arguments)]
pub(crate) async fn predict_and_emit(
    client: &PredictClient,
    fa: &clair3_connect_core::full_alignment::matrix::FaMatrix,
    candidates: &[u64],
    region_contig: &str,
    model_id: &str,
    matrix_depth: usize,
    batch_size: usize,
    output_config: &OutputWithConfig,
    predict_sem: &std::sync::Arc<tokio::sync::Semaphore>,
) -> Result<Vec<(u64, String)>> {
    let elements_per_row: usize = matrix_depth * fa.no_of_positions * fa.channel_size;
    let n = fa.candidate_num;
    const OUT_COLS: usize = 90;
    let mut y_all: Vec<f32> = vec![0.0; n * OUT_COLS];

    let mut handles = Vec::new();
    for row in (0..n).step_by(batch_size) {
        let end = (row + batch_size).min(n);
        let batch = end - row;
        let slice = &fa.matrix[row * elements_per_row..end * elements_per_row];
        let req = make_request_i8(
            &format!("fa_batch{}", row),
            model_id,
            vec![
                batch as i64,
                matrix_depth as i64,
                fa.no_of_positions as i64,
                fa.channel_size as i64,
            ],
            slice,
        );
        let mut c = client.clone();
        let sem = predict_sem.clone();
        handles.push(tokio::spawn(async move {
            let _permit = sem.acquire_owned().await.expect("predict semaphore closed");
            let resp = c
                .predict(req)
                .await
                .map_err(|e| anyhow!("fa predict batch {}: {}", row, e))?;
            Ok::<(usize, usize, Vec<f32>), anyhow::Error>((row, batch, decode_f32(&resp.y_payload)))
        }));
    }
    for h in handles {
        let (row, batch, y_batch) = h.await.map_err(|e| anyhow!("fa predict join: {e}"))??;
        let expected = batch * OUT_COLS;
        if y_batch.len() != expected {
            return Err(anyhow!(
                "fa Y batch {} unexpected size: got {}, want {}",
                row,
                y_batch.len(),
                expected
            ));
        }
        y_all[row * OUT_COLS..row * OUT_COLS + y_batch.len()].copy_from_slice(&y_batch);
    }

    let use_fa_path = output_config.add_indel_length;
    let mut emitted: Vec<(u64, String)> = Vec::with_capacity(n);
    for i in 0..n {
        let r = i * 90;
        let gt21 = &y_all[r..r + 21];
        let geno = &y_all[r + 21..r + 24];
        let vlen1 = &y_all[r + 24..r + 57];
        let vlen2 = &y_all[r + 57..r + 90];
        let candidate_1based = candidates[i] + 1;
        let alt_info = &fa.alt_info_strings[i];
        let mut parts = alt_info.splitn(4, '-');
        let _pos_field = parts.next().unwrap_or("0");
        let depth_field = parts.next().unwrap_or("0");
        let ref_byte = parts.next().and_then(|s| s.bytes().next()).unwrap_or(b'N');
        let payload_after = parts.next().unwrap_or("");
        let ref_char = ref_byte as char;
        let chr_pos_seq = format!("{}:{}:{}", region_contig, candidate_1based, ref_char);
        let alt_info_for_output = format!("{}-{}", depth_field, payload_after);
        let inputs_row = OutputWithInputs {
            chr_pos_seq: &chr_pos_seq,
            alt_info: &alt_info_for_output,
            gt21_probabilities: gt21,
            genotype_probabilities: geno,
            vlen_heads: if use_fa_path { Some((vlen1, vlen2)) } else { None },
        };
        match output_with(inputs_row, output_config) {
            Ok(Some(row_text)) => emitted.push((candidate_1based, row_text)),
            Ok(None) => {}
            Err(e) => tracing::warn!("fa output_with row {}: {e}", i),
        }
    }
    Ok(emitted)
}

pub(crate) async fn connect_predict_client(
    server_endpoint: &str,
    tls: Option<&clair3_connect_transport::tls::ClientTlsOptions>,
    aead_key: Option<&clair3_connect_transport::aead::AeadKey>,
) -> Result<PredictClient> {
    let mut client = match tls {
        Some(tls) => clair3_connect_transport::predict_client::PredictClient::connect_tls_with_retry(
            server_endpoint.to_string(),
            tls,
            clair3_connect_transport::predict_client::RetryPolicy::default(),
        )
        .await
        .context("connect_tls to server")?,
        None => PredictClient::connect(server_endpoint.to_string())
            .await
            .context("connect to server")?,
    };
    if let Some(key) = aead_key {
        client.set_aead_key(key.clone());
    }
    if let Some(secs) = std::env::var("CLAIR3_CONNECT_PREDICT_TIMEOUT_SECS")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
    {
        client.set_retry_policy(clair3_connect_transport::predict_client::RetryPolicy {
            per_attempt_timeout: Some(std::time::Duration::from_secs(secs)),
            ..clair3_connect_transport::predict_client::RetryPolicy::default()
        });
    }
    Ok(client)
}

pub async fn run_fa_pipeline<W: Write>(
    inputs: &FaPipelineInputs<'_>,
    out: &mut W,
) -> Result<()> {
    if inputs.candidates.is_empty() {

        let header = build_vcf_header(&VcfHeaderInputs {
            reference_file_path: Some(inputs.fasta_path),
            cmdline_file: None,
            sample_name: &inputs.sample_name,
            clair3_version: &inputs.clair3_version,
            gvcf: false,
        })?;
        out.write_all(header.text.as_bytes())?;
        out.write_all(b"\n")?;
        return Ok(());
    }

    let mut client = connect_predict_client(
        &inputs.server_endpoint,
        inputs.tls.as_ref(),
        inputs.aead_key.as_ref(),
    )
    .await?;

    let mut emitted: Vec<(u64, String)> = Vec::new();
    if inputs.low_mem {
        let window_bp = if inputs.fa_window_bp == 0 {
            2_000_000
        } else {
            inputs.fa_window_bp
        };

        const MARGIN: u64 = 64;
        let cands = &inputs.candidates;
        let mut i = 0usize;
        let mut n_windows = 0usize;
        while i < cands.len() {
            let win_start = cands[i];
            let mut j = i;
            while j < cands.len() && cands[j] < win_start.saturating_add(window_bp) {
                j += 1;
            }
            let win_cands = &cands[i..j];
            let region = Region {
                contig: inputs.region.contig.clone(),
                start: win_cands[0].saturating_sub(MARGIN),
                end: win_cands[win_cands.len() - 1] + MARGIN + 1,
            };
            let mut rows = fa_window_rows(&mut client, win_cands, &region, inputs).await?;
            emitted.append(&mut rows);
            n_windows += 1;
            i = j;
        }
        tracing::info!(
            "fa low-mem: {} candidates in {} windows of {} bp",
            cands.len(),
            n_windows,
            window_bp
        );
    } else {
        let mut rows =
            fa_window_rows(&mut client, &inputs.candidates, &inputs.region, inputs).await?;
        emitted.append(&mut rows);
    }
    emitted.sort_by_key(|(p, _)| *p);

    let header = build_vcf_header(&VcfHeaderInputs {
        reference_file_path: Some(inputs.fasta_path),
        cmdline_file: None,
        sample_name: &inputs.sample_name,
        clair3_version: &inputs.clair3_version,
        gvcf: false,
    })?;
    out.write_all(header.text.as_bytes())?;
    out.write_all(b"\n")?;
    for (_, r) in &emitted {
        out.write_all(r.as_bytes())?;
    }
    Ok(())
}
