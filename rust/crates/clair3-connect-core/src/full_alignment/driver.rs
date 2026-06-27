use super::candidates::FaCandidateContext;
use super::cigar_walk::{process_read, process_read_projected, PosAltInfo};
use super::constants::MIN_HAPLOTAG_MQ;
use super::haplotag::haplotag_read;
use super::matrix::{build_fa_matrix, build_fa_matrix_with_dwell, FaMatrix};
use super::phase_block_index::Variant;
use super::read_iter::{for_each_read, iterate_reads, iterate_reads_with_dwell};
use anyhow::{Context, Result};
use rust_htslib::faidx;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct FaInputs<'a> {

    pub bam_path: &'a Path,

    pub fasta_path: &'a Path,

    pub region: (String, i64, i64),

    pub candidates: &'a [u64],

    pub min_mq: u8,

    pub matrix_depth: usize,

    pub max_indel_length: usize,
}

pub fn calculate_full_alignment_no_haplotag(inputs: &FaInputs<'_>) -> Result<FaMatrix> {
    calculate_full_alignment(inputs, None, false)
}

pub fn calculate_full_alignment_streaming(inputs: &FaInputs<'_>) -> Result<FaMatrix> {
    let ctx = FaCandidateContext::build(inputs.candidates);
    let mut pos_alt_info: Vec<PosAltInfo> = vec![PosAltInfo::default(); ctx.candidate_num()];

    let (contig, region_start, region_end) =
        (inputs.region.0.as_str(), inputs.region.1, inputs.region.2);
    let expand: i64 = 2_000_000;
    let ref_start = (region_start - expand).max(0) as u64;
    let ref_end = (region_end + expand) as u64;
    let fasta = faidx::Reader::from_path(inputs.fasta_path).context("open fasta")?;
    let raw = fasta
        .fetch_seq(contig, ref_start as usize, ref_end as usize - 1)
        .context("fetch_seq for FA")?;
    let ref_seq: Vec<u8> = raw.iter().map(|b| b.to_ascii_uppercase()).collect();

    let mut projected = Vec::new();
    for_each_read(
        inputs.bam_path,
        &inputs.region,
        inputs.min_mq,
        false,
        |raw| {
            if let Some(pr) = process_read_projected(raw, &ctx, &mut pos_alt_info) {
                projected.push(pr);
            }
        },
    )
    .context("for_each_read (streaming FA)")?;

    Ok(build_fa_matrix(
        &ctx,
        &ref_seq,
        ref_start,
        &projected,
        &pos_alt_info,
        inputs.matrix_depth,
        inputs.max_indel_length,
    ))
}

pub fn calculate_full_alignment_with_haplotag(
    inputs: &FaInputs<'_>,
    variants: &[Variant],
) -> Result<FaMatrix> {
    calculate_full_alignment(inputs, Some(variants), false)
}

pub fn calculate_full_alignment_with_dwell(
    inputs: &FaInputs<'_>,
    variants: &[Variant],
) -> Result<FaMatrix> {
    calculate_full_alignment(inputs, Some(variants), true)
}

fn calculate_full_alignment(
    inputs: &FaInputs<'_>,
    variants: Option<&[Variant]>,
    enable_dwell_time: bool,
) -> Result<FaMatrix> {
    let ctx = FaCandidateContext::build(inputs.candidates);

    let reads = if enable_dwell_time {
        iterate_reads_with_dwell(inputs.bam_path, &inputs.region, inputs.min_mq)
    } else {
        iterate_reads(inputs.bam_path, &inputs.region, inputs.min_mq)
    }
    .context("iterate_reads")?;

    let mut pos_alt_info: Vec<PosAltInfo> = vec![PosAltInfo::default(); ctx.candidate_num()];

    let (contig, region_start, region_end) =
        (inputs.region.0.as_str(), inputs.region.1, inputs.region.2);
    let expand: i64 = 2_000_000;
    let ref_start = (region_start - expand).max(0) as u64;
    let ref_end = (region_end + expand) as u64;
    let fasta = faidx::Reader::from_path(inputs.fasta_path).context("open fasta")?;
    let raw = fasta
        .fetch_seq(contig, ref_start as usize, ref_end as usize - 1)
        .context("fetch_seq for FA")?;
    let ref_seq: Vec<u8> = raw.iter().map(|b| b.to_ascii_uppercase()).collect();

    let mut processed = Vec::with_capacity(reads.len());
    let mut variant_current_pos: usize = 0;
    for mut r in reads {

        if let Some(vs) = variants {

            while variant_current_pos < vs.len()
                && vs[variant_current_pos].position < r.read_start
            {
                variant_current_pos += 1;
            }

            if r.raw_mq >= MIN_HAPLOTAG_MQ {
                let (hap, new_read_end) = haplotag_read(
                    vs,
                    variant_current_pos,
                    r.read_start,
                    &r.cigar,
                    &r.seqi,
                    &ref_seq,
                    ref_start as i64,
                );
                r.haplotype = hap;
                r.read_end = new_read_end;
            }
        }

        if let Some(rwpi) = process_read(r, &ctx, &mut pos_alt_info) {
            processed.push(rwpi);
        }
    }

    Ok(if enable_dwell_time {
        build_fa_matrix_with_dwell(
            &ctx,
            &ref_seq,
            ref_start,
            &processed,
            &pos_alt_info,
            inputs.matrix_depth,
            inputs.max_indel_length,
        )
    } else {
        build_fa_matrix(
            &ctx,
            &ref_seq,
            ref_start,
            &processed,
            &pos_alt_info,
            inputs.matrix_depth,
            inputs.max_indel_length,
        )
    })
}

#[cfg(test)]
mod streaming_parity_tests {
    use super::*;
    use std::path::PathBuf;

    fn fixtures() -> Option<(PathBuf, PathBuf)> {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .to_path_buf();
        let bam = root.join("test_data/ont_chr20_demo/HG003_chr20_demo.bam");
        let fa = root.join("test_data/ont_chr20_demo/GRCh38_no_alt_chr20.fa");
        if bam.exists() && fa.exists() {
            Some((bam, fa))
        } else {
            None
        }
    }

    #[test]
    fn streaming_matches_batch_byte_for_byte() {
        let Some((bam, fa)) = fixtures() else {
            eprintln!("skip: chr20 demo fixture missing");
            return;
        };

        let candidates: Vec<u64> = (130_500..139_500).step_by(97).collect();
        let inputs = FaInputs {
            bam_path: &bam,
            fasta_path: &fa,
            region: ("chr20".to_string(), 130_000, 140_000),
            candidates: &candidates,
            min_mq: 5,
            matrix_depth: 89,
            max_indel_length: 50,
        };
        let batch = calculate_full_alignment_no_haplotag(&inputs).expect("batch build");
        let stream = calculate_full_alignment_streaming(&inputs).expect("streaming build");
        assert_eq!(batch.candidate_num, stream.candidate_num, "candidate_num");
        assert_eq!(batch.matrix_depth, stream.matrix_depth);
        assert_eq!(batch.no_of_positions, stream.no_of_positions);
        assert_eq!(batch.channel_size, stream.channel_size);
        assert_eq!(
            batch.matrix.len(),
            stream.matrix.len(),
            "matrix length differs"
        );
        assert!(
            batch.matrix == stream.matrix,
            "streaming matrix differs from batch (first diff at {:?})",
            batch
                .matrix
                .iter()
                .zip(&stream.matrix)
                .position(|(a, b)| a != b)
        );
        assert_eq!(
            batch.alt_info_strings, stream.alt_info_strings,
            "streaming alt_info differs from batch"
        );

        assert!(
            batch.matrix.iter().any(|&v| v != 0),
            "fixture slice had no coverage — test would be vacuous"
        );
    }
}
