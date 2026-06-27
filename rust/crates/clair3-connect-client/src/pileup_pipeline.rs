use anyhow::{anyhow, Context, Result};
use clair3_connect_core::hts::region::Region;
use clair3_connect_core::pileup::candidates::{window_candidates, WindowedCandidate};
use clair3_connect_core::pileup::chunk_stitch::split_on_gaps;
use clair3_connect_core::pileup::kernel::{calculate_pileup, PileupParams};
use clair3_connect_core::vcf::header::{build_vcf_header, VcfHeaderInputs};
use clair3_connect_core::vcf::output_with::{output_with, OutputWithConfig, OutputWithInputs};
use clair3_connect_transport::codec::decode_f32;
use clair3_connect_transport::predict_client::{make_request_i32, PredictClient};
use std::io::Write;
use std::path::Path;

pub fn expand_user_region(contig: &str, ctg_start_1based: u64, ctg_end_1based: u64) -> Region {
    const NO_OF_POSITIONS: u64 = 33;
    let extend_start_1b = ctg_start_1based.saturating_sub(NO_OF_POSITIONS).max(1);
    let extend_end_1b = ctg_end_1based + NO_OF_POSITIONS;

    let start_0b = extend_start_1b.saturating_sub(2);
    let end_0b = extend_end_1b;
    Region {
        contig: contig.to_string(),
        start: start_0b,
        end: end_0b,
    }
}

#[derive(Debug, Clone)]
pub struct PileupPipelineInputs<'a> {

    pub bam_path: &'a Path,

    pub fasta_path: &'a Path,

    pub region: Region,

    pub server_endpoint: String,

    pub model_id: String,

    pub sample_name: String,

    pub clair3_version: String,

    pub max_depth: i64,

    pub pileup_params: PileupParams,

    pub output_config: OutputWithConfig,

    pub batch_size: usize,

    pub aead_key: Option<clair3_connect_transport::aead::AeadKey>,

    pub tls: Option<clair3_connect_transport::tls::ClientTlsOptions>,
}

struct EmittedRow {
    pos: u64,
    row: String,
}

pub async fn run_pileup_pipeline<W: Write>(
    inputs: &PileupPipelineInputs<'_>,
    out: &mut W,
) -> Result<()> {

    let pileup = calculate_pileup(
        &inputs.region,
        inputs.bam_path,
        inputs.fasta_path,
        &inputs.pileup_params,
    )
    .context("calculate_pileup")?;
    tracing::info!(
        "calculate_pileup: {} chunks, {} candidates",
        pileup.chunks.len(),
        pileup.candidates.len()
    );

    let mut stitched = Vec::new();
    for chunk in pileup.chunks {
        let mut split = split_on_gaps(chunk);
        stitched.append(&mut split);
    }
    let windowed: Vec<WindowedCandidate> = window_candidates(&stitched, &pileup.candidates);
    tracing::info!("windowed: {} candidates survived", windowed.len());
    if windowed.is_empty() {

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

    let mut scaled_x: Vec<i32> = Vec::with_capacity(windowed.len() * 33 * 18);
    for w in &windowed {
        let alt_payload = &w.alt_info;
        let depth: i64 = alt_payload
            .split('-')
            .next()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let max_depth = inputs.max_depth;
        if depth > 0 && depth > max_depth * 3 / 2 {

            let scale: f64 = depth as f64 / max_depth as f64;
            for v in &w.tensor {

                let scaled = (*v as f64) / scale;
                scaled_x.push(scaled as i32);
            }
        } else {
            scaled_x.extend_from_slice(&w.tensor);
        }
    }

    let mut client = match inputs.tls.as_ref() {
        Some(tls) => clair3_connect_transport::predict_client::PredictClient::connect_tls_with_retry(
            inputs.server_endpoint.clone(),
            tls,
            clair3_connect_transport::predict_client::RetryPolicy::default(),
        )
        .await
        .context("connect_tls to server")?,
        None => PredictClient::connect(inputs.server_endpoint.clone())
            .await
            .context("connect to server")?,
    };
    if let Some(key) = inputs.aead_key.as_ref() {
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

    let elements_per_row: usize = 33 * 18;
    let n = windowed.len();
    const OUT_COLS: usize = 24;

    const PILEUP_PREDICT_CONCURRENCY: usize = 8;
    let mut y_all: Vec<f32> = vec![0.0; n * OUT_COLS];
    let batch_rows: Vec<usize> = (0..n).step_by(inputs.batch_size).collect();
    for wave in batch_rows.chunks(PILEUP_PREDICT_CONCURRENCY) {
        let mut handles = Vec::with_capacity(wave.len());
        for &row in wave {
            let end = (row + inputs.batch_size).min(n);
            let batch = end - row;
            let slice = &scaled_x[row * elements_per_row..end * elements_per_row];

            let req = make_request_i32(
                &format!("batch{}", row),
                &inputs.model_id,
                vec![batch as i64, 33, 18],
                slice,
            );
            let mut c = client.clone();
            handles.push(tokio::spawn(async move {
                let resp = c
                    .predict(req)
                    .await
                    .map_err(|e| anyhow!("predict batch {}: {}", row, e))?;
                Ok::<(usize, usize, Vec<f32>), anyhow::Error>((
                    row,
                    batch,
                    decode_f32(&resp.y_payload),
                ))
            }));
        }
        for h in handles {
            let (row, batch, y_batch) =
                h.await.map_err(|e| anyhow!("predict join: {e}"))??;
            let expected = batch * OUT_COLS;
            if y_batch.len() != expected {
                return Err(anyhow!(
                    "Y batch {} unexpected size: got {}, want {}",
                    row,
                    y_batch.len(),
                    expected
                ));
            }
            y_all[row * OUT_COLS..row * OUT_COLS + y_batch.len()].copy_from_slice(&y_batch);
        }
    }

    let mut emitted: Vec<EmittedRow> = Vec::with_capacity(windowed.len());
    for (i, w) in windowed.iter().enumerate() {
        let r = i * 24;
        let gt21 = &y_all[r..r + 21];
        let geno = &y_all[r + 21..r + 24];
        let chr_pos_seq = &w.pos_info;
        let alt_info_payload = &w.alt_info;
        let inputs_row = OutputWithInputs {
            chr_pos_seq,
            alt_info: alt_info_payload,
            gt21_probabilities: gt21,
            genotype_probabilities: geno,
            vlen_heads: None,
        };
        match output_with(inputs_row, &inputs.output_config) {
            Ok(Some(row_text)) => {

                let pos: u64 = row_text
                    .split('\t')
                    .nth(1)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
                emitted.push(EmittedRow {
                    pos,
                    row: row_text,
                });
            }
            Ok(None) => {}
            Err(e) => {
                tracing::warn!("output_with row {}: {e}", i);
            }
        }
    }

    emitted.sort_by_key(|e| e.pos);

    let header = build_vcf_header(&VcfHeaderInputs {
        reference_file_path: Some(inputs.fasta_path),
        cmdline_file: None,
        sample_name: &inputs.sample_name,
        clair3_version: &inputs.clair3_version,
        gvcf: false,
    })?;
    out.write_all(header.text.as_bytes())?;
    out.write_all(b"\n")?;
    for er in &emitted {
        out.write_all(er.row.as_bytes())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expand_user_region_matches_python_effective() {

        let r = expand_user_region("chr20", 100_000, 300_000);
        assert_eq!(r.start, 99_965);
        assert_eq!(r.end, 300_033);
    }

    #[test]
    fn expand_user_region_clamps_low() {
        let r = expand_user_region("chr1", 10, 100);

        assert_eq!(r.start, 0);
        assert_eq!(r.end, 133);
    }
}
