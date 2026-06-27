use crate::chunk_planner::ChunkRegion;
use crate::pileup_pipeline::{
    expand_user_region, run_pileup_pipeline, PileupPipelineInputs,
};
use anyhow::{anyhow, Context, Result};
use clair3_connect_core::manifest::{
    read_manifest_with_warnings, verify_manifest, ManifestEntry, ManifestWriter,
};
use clair3_connect_core::pileup::kernel::PileupParams;
use clair3_connect_core::vcf::header::{build_vcf_header, VcfHeaderInputs};
use clair3_connect_core::vcf::output_with::OutputWithConfig;
use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;
use tokio::sync::Semaphore;

#[derive(Debug, Clone)]
pub struct ParallelPileupConfig {

    pub bam_path: PathBuf,

    pub fasta_path: PathBuf,

    pub server_endpoint: String,

    pub model_id: String,

    pub sample_name: String,

    pub clair3_version: String,

    pub max_depth: i64,

    pub pileup_params: PileupParams,

    pub output_config: OutputWithConfig,

    pub batch_size: usize,

    pub max_concurrent: usize,

    pub manifest_dir: Option<PathBuf>,

    pub verify_resume: bool,

    pub aead_key: Option<clair3_connect_transport::aead::AeadKey>,

    pub tls: Option<clair3_connect_transport::tls::ClientTlsOptions>,
}

pub async fn run_parallel_pileup<W: Write>(
    chunks: &[ChunkRegion],
    config: &ParallelPileupConfig,
    out: &mut W,
) -> Result<ParallelStats> {
    if chunks.is_empty() {
        let header = build_vcf_header(&VcfHeaderInputs {
            reference_file_path: Some(&config.fasta_path),
            cmdline_file: None,
            sample_name: &config.sample_name,
            clair3_version: &config.clair3_version,
            gvcf: false,
        })?;
        out.write_all(header.text.as_bytes())?;
        out.write_all(b"\n")?;
        return Ok(ParallelStats::default());
    }

    let mut already_done: HashMap<String, ManifestEntry> = HashMap::new();
    let writer = if let Some(dir) = &config.manifest_dir {
        let (entries, warnings) =
            read_manifest_with_warnings(dir).context("read manifest")?;
        for w in &warnings {
            tracing::warn!("manifest: {:?}", w);
        }

        if config.verify_resume && !entries.is_empty() {
            verify_manifest(dir, &entries).context("verify resume manifest")?;
            tracing::info!(
                "resume: verified {} chunk vcfs against manifest sha256s",
                entries.len()
            );
        }
        for e in entries {
            already_done.insert(e.chunk_id.clone(), e);
        }
        tracing::info!(
            "resume: {} chunks already in manifest at {}",
            already_done.len(),
            dir.display()
        );
        Some(Arc::new(
            ManifestWriter::open(dir).context("open manifest writer")?,
        ))
    } else {
        None
    };

    let sem = Arc::new(Semaphore::new(config.max_concurrent.max(1)));
    let mut tasks: Vec<tokio::task::JoinHandle<Result<ChunkResult>>> = Vec::with_capacity(chunks.len());

    for chunk in chunks {

        let chunk_id = ManifestEntry::make_chunk_id(
            &chunk.contig,
            chunk.start_1based,
            chunk.end_1based,
        );
        if let Some(entry) = already_done.get(&chunk_id) {
            let dir = config.manifest_dir.as_ref().expect("manifest_dir set above");
            let path = dir.join(&entry.vcf_path);
            let text = std::fs::read_to_string(&path).context("read manifested chunk vcf")?;
            tracing::info!("resume: reusing chunk {} ({} bytes)", chunk_id, text.len());

            tasks.push(tokio::spawn(async move {
                Ok(ChunkResult {
                    text,
                    resumed: true,
                })
            }));
            continue;
        }
        let permit = sem.clone().acquire_owned().await.expect("sem acquire");
        let chunk = chunk.clone();
        let bam = config.bam_path.clone();
        let fa = config.fasta_path.clone();
        let server = config.server_endpoint.clone();
        let model = config.model_id.clone();
        let sample = config.sample_name.clone();
        let version = config.clair3_version.clone();
        let pileup_params = config.pileup_params.clone();
        let output_config = config.output_config;
        let max_depth = config.max_depth;
        let batch_size = config.batch_size;
        let aead_key = config.aead_key.clone();
        let tls = config.tls.clone();
        let writer_clone = writer.clone();

        tasks.push(tokio::spawn(async move {
            let _permit = permit;
            let mut buf: Vec<u8> = Vec::new();
            let inputs = PileupPipelineInputs {
                bam_path: &bam,
                fasta_path: &fa,
                region: expand_user_region(
                    &chunk.contig,
                    chunk.start_1based,
                    chunk.end_1based,
                ),
                server_endpoint: server,
                model_id: model,
                sample_name: sample,
                clair3_version: version,
                max_depth,
                pileup_params,
                output_config,
                batch_size,
                aead_key,
                tls,
            };
            run_pileup_pipeline(&inputs, &mut buf)
                .await
                .map_err(|e| anyhow!("chunk {}.{}: {}", chunk.contig, chunk.chunk_idx, e))?;
            let text = String::from_utf8(buf).map_err(|e| anyhow!("chunk utf8: {e}"))?;

            if let Some(w) = writer_clone {
                let vcf_path =
                    format!("chunk-{:05}-{}.vcf", chunk.chunk_idx, chunk.contig);
                let sha = w
                    .store_chunk_vcf(&vcf_path, &text)
                    .map_err(|e| anyhow!("store chunk vcf: {e}"))?;
                let row_count = text
                    .lines()
                    .filter(|l| !l.starts_with('#') && !l.is_empty())
                    .count() as u64;
                let entry = ManifestEntry {
                    chunk_id: ManifestEntry::make_chunk_id(
                        &chunk.contig,
                        chunk.start_1based,
                        chunk.end_1based,
                    ),
                    contig: chunk.contig.clone(),
                    start_1based: chunk.start_1based,
                    end_1based: chunk.end_1based,
                    row_count,
                    vcf_sha256: sha,
                    vcf_path,
                    source_phase: "pileup".to_string(),
                    completed_at_utc: chrono_utc_now(),
                };
                w.record_chunk(&entry)
                    .map_err(|e| anyhow!("record chunk: {e}"))?;
            }

            Ok(ChunkResult {
                text,
                resumed: false,
            })
        }));
    }

    let mut all_rows: Vec<(String, u64, String)> = Vec::new();
    let mut chunks_completed = 0_usize;
    let mut chunks_resumed = 0_usize;
    let mut chunks_failed = 0_usize;
    let mut first_error: Option<String> = None;
    for h in tasks {
        match h.await {
            Ok(Ok(r)) => {
                chunks_completed += 1;
                if r.resumed {
                    chunks_resumed += 1;
                }
                for line in r.text.lines() {
                    if line.starts_with('#') || line.is_empty() {
                        continue;
                    }
                    let cols: Vec<&str> = line.split('\t').collect();
                    if cols.len() < 2 {
                        continue;
                    }
                    let pos: u64 = match cols[1].parse() {
                        Ok(p) => p,
                        Err(_) => continue,
                    };
                    all_rows.push((cols[0].to_string(), pos, line.to_string()));
                }
            }
            Ok(Err(e)) => {
                chunks_failed += 1;
                if first_error.is_none() {
                    first_error = Some(e.to_string());
                }
            }
            Err(e) => {
                chunks_failed += 1;
                if first_error.is_none() {
                    first_error = Some(format!("join error: {e}"));
                }
            }
        }
    }
    if let Some(e) = first_error {
        return Err(anyhow!("{} chunks failed; first: {}", chunks_failed, e));
    }

    all_rows.sort_by(|a, b| {
        major_contig_rank(&a.0)
            .cmp(&major_contig_rank(&b.0))
            .then_with(|| a.0.cmp(&b.0))
            .then(a.1.cmp(&b.1))
    });

    let header = build_vcf_header(&VcfHeaderInputs {
        reference_file_path: Some(&config.fasta_path),
        cmdline_file: None,
        sample_name: &config.sample_name,
        clair3_version: &config.clair3_version,
        gvcf: false,
    })?;
    out.write_all(header.text.as_bytes())?;
    out.write_all(b"\n")?;
    for (_, _, row) in &all_rows {
        out.write_all(row.as_bytes())?;
        out.write_all(b"\n")?;
    }

    Ok(ParallelStats {
        chunks_total: chunks.len(),
        chunks_completed,
        chunks_resumed,
        chunks_failed,
        rows_emitted: all_rows.len(),
    })
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct ParallelStats {

    pub chunks_total: usize,

    pub chunks_completed: usize,

    pub chunks_resumed: usize,

    pub chunks_failed: usize,

    pub rows_emitted: usize,
}

struct ChunkResult {
    text: String,

    resumed: bool,
}

fn chrono_utc_now() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);

    let (year, month, day, hour, minute, second) = unix_to_ymdhms(secs);
    format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}Z",
        year, month, day, hour, minute, second
    )
}

fn unix_to_ymdhms(secs: i64) -> (i32, u32, u32, u32, u32, u32) {
    let secs = secs.max(0);
    let day = secs / 86_400;
    let rem = secs % 86_400;
    let hour = (rem / 3_600) as u32;
    let minute = ((rem % 3_600) / 60) as u32;
    let second = (rem % 60) as u32;

    let mut year = 1970_i32;
    let mut days_remaining = day;
    loop {
        let leap = is_leap(year);
        let yd = if leap { 366 } else { 365 };
        if days_remaining < yd {
            break;
        }
        days_remaining -= yd;
        year += 1;
    }
    let leap = is_leap(year);
    let month_days = [31, if leap { 29 } else { 28 }, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    let mut month = 1_u32;
    for &md in &month_days {
        if days_remaining < md {
            break;
        }
        days_remaining -= md;
        month += 1;
    }
    let day = (days_remaining + 1) as u32;
    (year, month, day, hour, minute, second)
}

fn is_leap(y: i32) -> bool {
    (y % 4 == 0 && y % 100 != 0) || y % 400 == 0
}

fn major_contig_rank(contig: &str) -> u32 {
    let stripped = contig.strip_prefix("chr").unwrap_or(contig);
    if let Ok(n) = stripped.parse::<u32>() {
        if (1..=22).contains(&n) {
            return n;
        }
    }
    match stripped {
        "X" => 23,
        "Y" => 24,
        "M" | "MT" => 25,
        _ => 100,
    }
}

#[cfg(test)]
mod tests {
    use super::major_contig_rank;

    #[test]
    fn major_contig_ranks_chr_numeric() {
        assert_eq!(major_contig_rank("chr1"), 1);
        assert_eq!(major_contig_rank("chr2"), 2);
        assert_eq!(major_contig_rank("chr10"), 10);
        assert_eq!(major_contig_rank("chr22"), 22);

        assert_eq!(major_contig_rank("1"), 1);
        assert_eq!(major_contig_rank("22"), 22);
    }

    #[test]
    fn major_contig_ranks_xy_m() {
        assert_eq!(major_contig_rank("chrX"), 23);
        assert_eq!(major_contig_rank("X"), 23);
        assert_eq!(major_contig_rank("chrY"), 24);
        assert_eq!(major_contig_rank("chrM"), 25);
        assert_eq!(major_contig_rank("chrMT"), 25);
    }

    #[test]
    fn major_contig_ranks_unknown_last() {
        assert_eq!(major_contig_rank("GL000226.1"), 100);
        assert_eq!(major_contig_rank("hs37d5"), 100);
    }

    #[test]
    fn major_contig_ordering_chr1_before_chr10_before_chr2() {

        let mut v = ["chr10", "chr2", "chr1", "chr20"];
        v.sort_by_key(|c| major_contig_rank(c));
        assert_eq!(v, ["chr1", "chr2", "chr10", "chr20"]);
    }
}
