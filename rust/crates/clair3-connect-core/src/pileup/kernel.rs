use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{Context, Result};
use rust_htslib::bam::Read;

use super::tables::*;
use crate::hts::bam_set::{BamSet, BAM_PILEUP_FILTER_FLAGS};
use crate::hts::fasta::fetch_window;
use crate::hts::region::Region;

#[derive(Debug, Clone)]
pub struct PileupParams {

    pub min_depth: u64,

    pub min_snp_af: f32,

    pub min_indel_af: f32,

    pub min_mapq: u8,

    pub max_indel_length: usize,

    pub call_snp_only: bool,

    pub max_depth: usize,

    pub gvcf: bool,

    pub call_ht: bool,
}

impl Default for PileupParams {
    fn default() -> Self {
        Self {
            min_depth: 2,
            min_snp_af: 0.08,
            min_indel_af: 0.15,
            min_mapq: 5,
            max_indel_length: 50,
            call_snp_only: false,
            max_depth: 144,
            gvcf: false,
            call_ht: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct PileupChunk {

    pub matrix: Vec<i64>,

    pub major: Vec<i64>,

    pub minor: Vec<i64>,

    pub n_cols: usize,
}

#[derive(Debug, Clone)]
pub struct AltInfo {

    pub pos_1based: i64,

    pub pos_info: String,

    pub alt_info: String,
}

#[derive(Debug)]
pub struct PileupResult {

    pub chunks: Vec<PileupChunk>,

    pub candidates: Vec<AltInfo>,

    pub gvcf_pos_ref_count: Vec<i64>,

    pub gvcf_pos_total_count: Vec<i64>,
}

pub fn calculate_pileup(
    region: &Region,
    bam_path: &Path,
    fasta_path: &Path,
    params: &PileupParams,
) -> Result<PileupResult> {
    let mut bam_set = BamSet::open(bam_path, fasta_path)?;

    let bgzf_threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
        .min(4);
    if bgzf_threads > 1 {
        let _ = bam_set.reader.set_threads(bgzf_threads);
    }
    let tid = bam_set.tid(&region.contig)?;
    let ref_window = fetch_window(&bam_set.fasta, &region.contig, region.start, region.end)?;

    bam_set
        .reader
        .fetch((tid, region.start as i64, region.end as i64))
        .with_context(|| {
            format!(
                "fetch BAM {}:{}-{}",
                region.contig, region.start, region.end
            )
        })?;

    let mut pileups = bam_set.reader.pileup();

    pileups.set_max_depth(1 << 20);

    let mut matrix: Vec<i64> = Vec::with_capacity((region.end - region.start) as usize * FEAT_LEN);
    let mut major: Vec<i64> = Vec::with_capacity((region.end - region.start) as usize);
    let mut minor: Vec<i64> = Vec::with_capacity((region.end - region.start) as usize);
    let mut candidates: Vec<AltInfo> = Vec::new();

    let mut pre_pos: u64 = 0;
    let mut contiguous_flanking_num: u64 = 0;

    let mut is_del_ins_encountered: u64 = 0;

    for plp in pileups {
        let pileup = plp.context("pileup iterator")?;
        if pileup.tid() != tid {
            continue;
        }
        let pos = pileup.pos() as u64;

        if pos < region.start {
            continue;
        }
        if pos >= region.end {
            break;
        }

        if pre_pos + 1 != pos || pre_pos == 0 {
            contiguous_flanking_num = 0;
        } else {
            contiguous_flanking_num += 1;
        }
        pre_pos = pos;

        let mut row = [0i64; FEAT_LEN];

        let mut dels_f: Vec<u64> = Vec::new();
        let mut dels_r: Vec<u64> = Vec::new();

        let mut ins_f: BTreeMap<Vec<u8>, u64> = BTreeMap::new();
        let mut ins_r: BTreeMap<Vec<u8>, u64> = BTreeMap::new();
        let mut ins_all: BTreeMap<Vec<u8>, u64> = BTreeMap::new();

        let mut depth: u64 = 0;

        for align in pileup.alignments() {
            let record = align.record();

            if (record.flags() & BAM_PILEUP_FILTER_FLAGS) != 0 {
                continue;
            }
            if record.mapq() < params.min_mapq {
                continue;
            }

            if align.is_refskip() {
                continue;
            }

            if let rust_htslib::bam::pileup::Indel::Del(d) = align.indel() {
                let d_usize = d as usize;
                if dels_f.len() < d_usize {
                    dels_f.resize(d_usize, 0);
                }
                if dels_r.len() < d_usize {
                    dels_r.resize(d_usize, 0);
                }
                if record.is_reverse() {
                    dels_r[d_usize - 1] += 1;
                } else {
                    dels_f[d_usize - 1] += 1;
                }
            }

            if align.is_del() {
                let col = if record.is_reverse() {
                    C3_REV_DEL
                } else {
                    C3_FWD_DEL
                };
                row[col] += 1;
                depth += 1;
            } else {

                let qpos = match align.qpos() {
                    Some(q) => q,
                    None => continue,
                };
                let seq = record.seq();
                let base_j = seq.encoded_base(qpos) as usize;
                let lookup_idx = if record.is_reverse() {
                    base_j + 16
                } else {
                    base_j
                };
                let base_i = NUM_TO_COUNT_BASE_CLAIR3[lookup_idx];

                if base_i < 0 {

                    if !matrix.is_empty() {
                        let last_idx = matrix.len() - 1;
                        matrix[last_idx] = matrix[last_idx].saturating_add(1);
                    }

                } else {
                    row[base_i as usize] += 1;
                }
                depth += 1;
            }

            if let rust_htslib::bam::pileup::Indel::Ins(d) = align.indel() {
                let d = d as usize;
                let qpos_start: Option<usize> = match align.qpos() {
                    Some(q) => Some(q + 1),
                    None => {

                        is_del_ins_encountered += 1;
                        match record.cigar().read_pos(pileup.pos(), false, true) {
                            Ok(Some(rp)) => Some(rp as usize),
                            _ => None,
                        }
                    }
                };
                if let Some(qstart) = qpos_start {
                    let seq = record.seq();
                    let mut ins_str: Vec<u8> = Vec::with_capacity(d);
                    for j in 0..d {
                        let nibble = seq.encoded_base(qstart + j) as usize;
                        ins_str.push(SEQ_NT16_STR[nibble]);
                    }
                    if record.is_reverse() {
                        *ins_r.entry(ins_str.clone()).or_insert(0) += 1;
                    } else {
                        *ins_f.entry(ins_str.clone()).or_insert(0) += 1;
                    }
                    *ins_all.entry(ins_str).or_insert(0) += 1;
                }
            }
        }

        let (fwd_del_all, fwd_del_best) = sum_max_vec(&dels_f);
        let (rev_del_all, rev_del_best) = sum_max_vec(&dels_r);
        row[C3_FWD_DEL_ALL] = fwd_del_all as i64;
        row[C3_FWD_DEL_BEST] = fwd_del_best as i64;
        row[C3_REV_DEL_ALL] = rev_del_all as i64;
        row[C3_REV_DEL_BEST] = rev_del_best as i64;
        let del_count = fwd_del_all + rev_del_all;

        let (fwd_ins_all, fwd_ins_best) = sum_max_map(&ins_f);
        let (rev_ins_all, rev_ins_best) = sum_max_map(&ins_r);
        row[C3_FWD_INS_ALL] = fwd_ins_all as i64;
        row[C3_FWD_INS_BEST] = fwd_ins_best as i64;
        row[C3_REV_INS_ALL] = rev_ins_all as i64;
        row[C3_REV_INS_BEST] = rev_ins_best as i64;
        let ins_count = fwd_ins_all + rev_ins_all;

        let ref_base = ref_window.base_at(pos).unwrap_or(b'N');
        let ref_offset_forward: usize = if ref_base.is_ascii_uppercase() {

            BASE2INDEX[(ref_base - b'A') as usize] as usize
        } else {
            0
        };
        let ref_offset_reverse = ref_offset_forward + REVERSE_POS_START;

        let mut forward_sum: i64 = 0;
        let mut reverse_sum: i64 = 0;
        let mut ref_count: i64 = 0;
        let mut alt_count: i64 = 0;
        let mut major_alt_base: u8 = 0;
        for i in 0..4 {
            forward_sum += row[i];
            reverse_sum += row[i + REVERSE_POS_START];
            if i == ref_offset_forward {
                ref_count = row[i] + row[i + REVERSE_POS_START];
            } else {
                let current_count = row[i] + row[i + REVERSE_POS_START];
                if current_count > alt_count {
                    alt_count = current_count;
                    major_alt_base = PLP_BASES_CLAIR3[i];
                }
            }
        }

        row[ref_offset_forward] = -forward_sum;
        row[ref_offset_reverse] = -reverse_sum;

        let depth_eff = depth.max(1) as i64;
        let pass_min_depth = depth_eff >= params.min_depth as i64;
        let pass_ref_base_in_acgt = matches!(ref_base, b'A' | b'C' | b'G' | b'T');
        let non_ref_base_majority =
            ref_count < alt_count || ref_count < ins_count as i64 || ref_count < del_count as i64;
        let ref_alt_equal_majority =
            ref_count > 0 && ref_count == alt_count && (ref_base as i32) < (major_alt_base as i32);
        let af_snp = alt_count as f32 / depth_eff as f32;
        let af_del = del_count as f32 / depth_eff as f32;
        let af_ins = ins_count as f32 / depth_eff as f32;
        let mut pass_af = if params.call_snp_only {
            af_snp >= params.min_snp_af
        } else {
            non_ref_base_majority
                || ref_alt_equal_majority
                || af_snp >= params.min_snp_af
                || af_del >= params.min_indel_af
                || af_ins >= params.min_indel_af
        };
        pass_af = pass_af && pass_min_depth && pass_ref_base_in_acgt;
        if !params.call_ht {
            pass_af = pass_af && contiguous_flanking_num >= FLANKING_BASE_NUM;
        }

        if pass_af {

            let mut alt_info = format!("{}-{}-{}-", pos + 1, depth_eff as u64, ref_base as char);

            for i in 0..4 {
                if i == ref_offset_forward {
                    continue;
                }
                let alt_sum = row[i] + row[i + REVERSE_POS_START];
                if alt_sum > 0 {
                    alt_info.push_str(&format!("X{} {} ", PLP_BASES_CLAIR3[i] as char, alt_sum));
                }
            }

            let offset = (pos - ref_window.start) as usize;
            let mut ref_depth: i64 = ref_count;
            for (i, &df) in dels_f.iter().enumerate() {
                let dr = dels_r.get(i).copied().unwrap_or(0);
                let d = df + dr;
                ref_depth -= d as i64;
                if d > 0 && (i + 1) <= params.max_indel_length {
                    let del_len = i + 1;
                    let seq_start = offset + 1;
                    let seq_end = seq_start + del_len;
                    if let Some(seq) = ref_window.seq.get(seq_start..seq_end) {
                        let seq_str = std::str::from_utf8(seq).unwrap_or("?");
                        alt_info.push_str(&format!("D{} {} ", seq_str, d));
                    } else {
                        tracing::warn!(
                            contig = %region.contig, pos, del_len,
                            ref_window_len = ref_window.seq.len(),
                            "ref window too short for D entry; alt_info will diverge from libclair3.so output"
                        );
                    }
                }
            }

            for (i, &dr) in dels_r.iter().enumerate().skip(dels_f.len()) {
                ref_depth -= dr as i64;
                if dr > 0 && (i + 1) <= params.max_indel_length {
                    let del_len = i + 1;
                    let seq_start = offset + 1;
                    let seq_end = seq_start + del_len;
                    if let Some(seq) = ref_window.seq.get(seq_start..seq_end) {
                        let seq_str = std::str::from_utf8(seq).unwrap_or("?");
                        alt_info.push_str(&format!("D{} {} ", seq_str, dr));
                    } else {
                        tracing::warn!(
                            contig = %region.contig, pos, del_len,
                            ref_window_len = ref_window.seq.len(),
                            "ref window too short for D entry (rev-only); alt_info will diverge from libclair3.so output"
                        );
                    }
                }
            }

            for (key, &val) in &ins_all {
                ref_depth -= val as i64;
                if key.len() <= params.max_indel_length {
                    let key_str = std::str::from_utf8(key).unwrap_or("?");
                    alt_info.push_str(&format!("I{}{} {} ", ref_base as char, key_str, val));
                }
            }

            if ref_depth > 0 {
                alt_info.push_str(&format!("R{} {} ", ref_base as char, ref_depth));
            }

            let pos_info = format!("{}:{}:{}", region.contig, pos + 1, ref_base as char);

            let depth_str = depth_eff.to_string();
            let payload_offset = format!("{}-{}-{}-", pos + 1, depth_str, ref_base as char).len();

            let payload = alt_info[payload_offset..].trim_end();
            let alt_info_field = format!("{}-{}", depth_str, payload);
            candidates.push(AltInfo {
                pos_1based: (pos + 1) as i64,
                pos_info,
                alt_info: alt_info_field,
            });
        }

        major.push(pos as i64);
        minor.push(0);
        matrix.extend_from_slice(&row);
    }

    if is_del_ins_encountered > 0 {
        tracing::debug!(
            count = is_del_ins_encountered,
            "is_del + Ins reads recovered via CIGAR-walk in {}:{}-{}",
            region.contig,
            region.start,
            region.end
        );
    }

    let n_cols = major.len();
    let single = PileupChunk {
        matrix,
        major,
        minor,
        n_cols,
    };

    let chunks = super::chunk_stitch::split_on_gaps(single);

    Ok(PileupResult {
        chunks,
        candidates,
        gvcf_pos_ref_count: Vec::new(),
        gvcf_pos_total_count: Vec::new(),
    })
}

fn sum_max_vec(v: &[u64]) -> (u64, u64) {
    let mut sum = 0u64;
    let mut best = 0u64;
    for &x in v {
        sum = sum.saturating_add(x);
        if x > best {
            best = x;
        }
    }
    (sum, best)
}

fn sum_max_map(m: &BTreeMap<Vec<u8>, u64>) -> (u64, u64) {
    let mut sum = 0u64;
    let mut best = 0u64;
    for &v in m.values() {
        sum = sum.saturating_add(v);
        if v > best {
            best = v;
        }
    }
    (sum, best)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};
    use std::collections::BTreeMap;

    #[test]
    fn sum_max_vec_basics() {
        assert_eq!(sum_max_vec(&[]), (0, 0));
        assert_eq!(sum_max_vec(&[3, 1, 4, 1, 5, 9, 2, 6]), (31, 9));
        assert_eq!(sum_max_vec(&[0, 0, 0]), (0, 0));
        assert_eq!(sum_max_vec(&[7]), (7, 7));
    }

    #[test]
    fn sum_max_map_basics() {
        let mut m: BTreeMap<Vec<u8>, u64> = BTreeMap::new();
        assert_eq!(sum_max_map(&m), (0, 0));
        m.insert(b"A".to_vec(), 3);
        m.insert(b"AC".to_vec(), 5);
        m.insert(b"AG".to_vec(), 1);
        assert_eq!(sum_max_map(&m), (9, 5));
    }

    #[test]
    fn is_del_ins_cigar_walk_offset() {
        let cigar = CigarString(vec![
            Cigar::Match(5),
            Cigar::Del(2),
            Cigar::Ins(3),
            Cigar::Match(3),
        ]);
        let view = cigar.into_view(100);
        let rp = view
            .read_pos(105, false, true)
            .expect("read_pos OK")
            .expect("Some(rp)");

        assert_eq!(rp, 5, "expected insertion-start read pos");

        let qstart = rp as usize;
        let qend_excl = qstart + 3;
        assert_eq!((qstart, qend_excl), (5, 8));
    }
}
