use super::candidates::FaCandidateContext;
use super::cigar_walk::{PosAltInfo, PosInfo, ProjectedRead, ReadWithPosInfo};
use super::constants::{normalize_af, normalize_hap, FLANKING_BASE_NUM, NO_OF_POSITIONS};
use super::tables::{acgt_num_for_char, fa_feature_for_char};

pub trait FaReadView {

    fn read_start(&self) -> i64;

    fn read_end(&self) -> i64;

    fn haplotype(&self) -> u8;

    fn strand(&self) -> i8;

    fn mq(&self) -> i8;

    fn flanking_start(&self) -> usize;

    fn overlap_candidates_num(&self) -> usize;

    fn pos_info(&self) -> &[PosInfo];
}

impl FaReadView for ReadWithPosInfo {
    #[inline]
    fn read_start(&self) -> i64 {
        self.raw.read_start
    }
    #[inline]
    fn read_end(&self) -> i64 {
        self.raw.read_end
    }
    #[inline]
    fn haplotype(&self) -> u8 {
        self.raw.haplotype
    }
    #[inline]
    fn strand(&self) -> i8 {
        self.raw.strand
    }
    #[inline]
    fn mq(&self) -> i8 {
        self.raw.mq
    }
    #[inline]
    fn flanking_start(&self) -> usize {
        self.flanking_start
    }
    #[inline]
    fn overlap_candidates_num(&self) -> usize {
        self.overlap_candidates_num
    }
    #[inline]
    fn pos_info(&self) -> &[PosInfo] {
        &self.pos_info
    }
}

impl FaReadView for ProjectedRead {
    #[inline]
    fn read_start(&self) -> i64 {
        self.read_start
    }
    #[inline]
    fn read_end(&self) -> i64 {
        self.read_end
    }
    #[inline]
    fn haplotype(&self) -> u8 {
        self.haplotype
    }
    #[inline]
    fn strand(&self) -> i8 {
        self.strand
    }
    #[inline]
    fn mq(&self) -> i8 {
        self.mq
    }
    #[inline]
    fn flanking_start(&self) -> usize {
        self.flanking_start
    }
    #[inline]
    fn overlap_candidates_num(&self) -> usize {
        self.overlap_candidates_num
    }
    #[inline]
    fn pos_info(&self) -> &[PosInfo] {
        &self.pos_info
    }
}

#[derive(Debug)]
pub struct FaMatrix {

    pub matrix: Vec<i8>,

    pub alt_info_strings: Vec<String>,

    pub candidate_num: usize,

    pub matrix_depth: usize,

    pub no_of_positions: usize,

    pub channel_size: usize,
}

pub fn build_fa_matrix<R: FaReadView>(
    ctx: &FaCandidateContext,
    ref_seq: &[u8],
    ref_start: u64,
    reads: &[R],
    pos_alt_info: &[PosAltInfo],
    matrix_depth: usize,
    max_indel_length: usize,
) -> FaMatrix {
    build_fa_matrix_inner(ctx, ref_seq, ref_start, reads, pos_alt_info, matrix_depth, max_indel_length, false)
}

pub fn build_fa_matrix_with_dwell<R: FaReadView>(
    ctx: &FaCandidateContext,
    ref_seq: &[u8],
    ref_start: u64,
    reads: &[R],
    pos_alt_info: &[PosAltInfo],
    matrix_depth: usize,
    max_indel_length: usize,
) -> FaMatrix {
    build_fa_matrix_inner(ctx, ref_seq, ref_start, reads, pos_alt_info, matrix_depth, max_indel_length, true)
}

#[allow(clippy::too_many_arguments)]
fn build_fa_matrix_inner<R: FaReadView>(
    ctx: &FaCandidateContext,
    ref_seq: &[u8],
    ref_start: u64,
    reads: &[R],
    pos_alt_info: &[PosAltInfo],
    matrix_depth: usize,
    max_indel_length: usize,
    enable_dwell_time: bool,
) -> FaMatrix {
    let channel_size = if enable_dwell_time { 9 } else { 8 };
    let signal_channel_offset = channel_size - 1;
    let n_candidates = ctx.candidate_num();
    let offset_can = NO_OF_POSITIONS * matrix_depth * channel_size;
    let offset_row = NO_OF_POSITIONS * channel_size;
    let offset_col = channel_size;
    let mut matrix: Vec<i8> = vec![0; n_candidates * offset_can];
    let mut alt_info_strings: Vec<String> = Vec::with_capacity(n_candidates);

    #[derive(Default, Clone)]
    struct RowAltInfo {
        alt_base: u8,
        ins_bases: Vec<u8>,
        del_length: u32,
        has_alt_info: bool,
    }

    for (i, pai) in pos_alt_info.iter().enumerate().take(n_candidates) {
        let candidate = ctx.candidates[i];

        let start_pos_signed = candidate as i64 - FLANKING_BASE_NUM as i64;
        let start_pos: u64 = start_pos_signed.max(0) as u64;
        let end_pos: u64 = candidate + FLANKING_BASE_NUM as u64 + 1;

        let mut overlap: Vec<(usize, u8)> = Vec::new();
        for (j, r) in reads.iter().enumerate() {
            if r.read_start() as u64 >= end_pos {
                break;
            }
            if (r.read_end() as u64) <= start_pos {
                continue;
            }
            overlap.push((j, r.haplotype()));
        }

        if overlap.len() > matrix_depth {

            fn mix64(mut z: u64) -> u64 {
                z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^ (z >> 31)
            }
            overlap.sort_by_key(|&(ri, _)| mix64(ri as u64));
            overlap.truncate(matrix_depth);
        }
        overlap.sort_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));

        let read_num = overlap.len();
        let padding = matrix_depth - read_num;
        let prefix_pad = padding / 2;
        let suffix_pad = padding - prefix_pad;
        let mut row_for_index: Vec<Option<usize>> = Vec::with_capacity(matrix_depth);
        for _ in 0..prefix_pad {
            row_for_index.push(None);
        }
        for (ri, _) in &overlap {
            row_for_index.push(Some(*ri));
        }
        for _ in 0..suffix_pad {
            row_for_index.push(None);
        }
        debug_assert_eq!(row_for_index.len(), matrix_depth);

        let mut alt_info_rows: Vec<RowAltInfo> = vec![RowAltInfo::default(); matrix_depth];

        for (d, rfi) in row_for_index.iter().enumerate() {
            let Some(read_index) = rfi else {
                continue;
            };
            let read = &reads[*read_index];
            let hap_v = normalize_hap(read.haplotype());
            let strand_v = read.strand();
            let mq_v = read.mq();

            for p in 0..NO_OF_POSITIONS {
                let cp = start_pos + p as u64;
                let is_center_pos = p == FLANKING_BASE_NUM;

                let Some(flanking_index) = ctx.flanking_idx(cp) else {
                    continue;
                };
                if flanking_index < read.flanking_start() {
                    continue;
                }
                let offset = flanking_index - read.flanking_start();
                if offset >= read.overlap_candidates_num() {
                    continue;
                }
                let pi = &read.pos_info()[offset];

                if pi.alt_base < 0 {
                    continue;
                }

                let ref_local = (cp as i64 - ref_start as i64) as usize;
                if ref_local >= ref_seq.len() {
                    continue;
                }
                let ref_base = ref_seq[ref_local].to_ascii_uppercase();
                let ref_v = fa_feature_for_char(ref_base);
                let bq_v = pi.bq;

                let _is_center_pos = is_center_pos;
                let alt_base_char = super::tables::seqi_to_char(pi.alt_base as u8);

                let mut alt_v = 0_i8;
                if !pi.ins_bases.is_empty() {
                    let ins_length = pi.ins_bases.len();
                    let remaining_cols = NO_OF_POSITIONS - p;
                    let max_ins_length = ins_length.min(remaining_cols);
                    for ins_idx in 0..ins_length {
                        let ins_alt_base = pi.ins_bases[ins_idx];
                        if ins_idx < max_ins_length && p < NO_OF_POSITIONS - 1 {
                            let ins_v = fa_feature_for_char(ins_alt_base);
                            let idx = i * offset_can
                                + d * offset_row
                                + (ins_idx + p) * offset_col
                                + 6;
                            matrix[idx] = ins_v;
                        }
                    }
                    if is_center_pos {
                        alt_info_rows[d].alt_base = alt_base_char;
                        alt_info_rows[d].ins_bases = pi.ins_bases.clone();
                        alt_info_rows[d].has_alt_info = true;
                    }
                    alt_v = fa_feature_for_char(b'I');
                } else if pi.del_length > 0 {
                    if is_center_pos {
                        alt_info_rows[d].del_length = pi.del_length;
                        alt_info_rows[d].has_alt_info = true;
                    }
                    alt_v = fa_feature_for_char(b'D');
                } else if ref_base != alt_base_char {
                    if is_center_pos {
                        alt_info_rows[d].alt_base = alt_base_char;
                        alt_info_rows[d].has_alt_info = true;
                    }
                    alt_v = fa_feature_for_char(alt_base_char);
                }

                let base = i * offset_can + d * offset_row + p * offset_col;
                matrix[base] = ref_v;
                matrix[base + 1] = alt_v;
                matrix[base + 2] = strand_v;
                matrix[base + 3] = mq_v;
                matrix[base + 4] = bq_v;
                matrix[base + 7] = hap_v;
                if enable_dwell_time {

                    let raw_signal = pi.signal_length;
                    matrix[base + signal_channel_offset] = raw_signal as i8;
                }
            }
        }

        let candidate_depth = pai.depth;
        for (d, row) in alt_info_rows.iter().enumerate().take(matrix_depth) {
            if !row.has_alt_info {
                continue;
            }
            let af_v: i8 = if !row.ins_bases.is_empty() {
                let count = *pos_alt_info[i]
                    .ins_counter
                    .get(&row.ins_bases)
                    .unwrap_or(&0);
                if count > 0 && candidate_depth > 0 {
                    normalize_af(count as f64 / candidate_depth as f64)
                } else {
                    0
                }
            } else if row.del_length > 0 {
                let count = *pos_alt_info[i]
                    .del_counter
                    .get(&row.del_length)
                    .unwrap_or(&0);
                if count > 0 && candidate_depth > 0 {
                    normalize_af(count as f64 / candidate_depth as f64)
                } else {
                    0
                }
            } else if row.alt_base != 0 {
                let bucket = acgt_num_for_char(row.alt_base);
                if (0..4).contains(&bucket) && candidate_depth > 0 {
                    let count = pai.acgt_count[bucket as usize];
                    normalize_af(count as f64 / candidate_depth as f64)
                } else {
                    0
                }
            } else {
                0
            };

            if af_v > 0 {
                for p in 0..NO_OF_POSITIONS {
                    let base = i * offset_can + d * offset_row + p * offset_col;
                    if matrix[base] != 0 {
                        matrix[base + 5] = af_v;
                    }
                }
            }
        }

        let center_ref_byte = {
            let local = (candidate as i64 - ref_start as i64) as usize;
            ref_seq.get(local).copied().unwrap_or(b'N').to_ascii_uppercase()
        };
        let mut ref_count: i64 = {
            let bucket = acgt_num_for_char(center_ref_byte);
            if (0..4).contains(&bucket) {
                pai.acgt_count[bucket as usize] as i64
            } else {
                0
            }
        };

        let mut s = String::new();
        s.push_str(&format!(
            "{}-{}-{}-",
            candidate + 1,
            candidate_depth,
            center_ref_byte as char
        ));

        const ACGT_LETTERS: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let ref_bucket = acgt_num_for_char(center_ref_byte);
        for (j, base) in ACGT_LETTERS.iter().enumerate() {
            if (j as i8) != ref_bucket && pai.acgt_count[j] > 0 {
                s.push_str(&format!(
                    "X{} {} ",
                    *base as char,
                    pai.acgt_count[j]
                ));
            }
        }

        for (key, val) in &pai.ins_counter {
            ref_count -= *val as i64;
            if key.len() <= max_indel_length {
                s.push_str(&format!(
                    "I{}{} {} ",
                    center_ref_byte as char,
                    std::str::from_utf8(key).unwrap_or("?"),
                    val
                ));
            }
        }

        for (key, val) in &pai.del_counter {
            ref_count -= *val as i64;
            if (*key as usize) <= max_indel_length {

                let key_len = *key as usize;
                let local_start = (candidate as i64 - ref_start as i64) as usize + 1;
                let local_end = local_start + key_len;
                if local_end <= ref_seq.len() {
                    let deleted_bases = std::str::from_utf8(&ref_seq[local_start..local_end])
                        .unwrap_or("?")
                        .to_ascii_uppercase();
                    s.push_str(&format!("D{} {} ", deleted_bases, val));
                }
            }
        }
        if ref_count > 0 {
            s.push_str(&format!("R{} {} ", center_ref_byte as char, ref_count));
        }
        alt_info_strings.push(s);
    }

    FaMatrix {
        matrix,
        alt_info_strings,
        candidate_num: n_candidates,
        matrix_depth,
        no_of_positions: NO_OF_POSITIONS,
        channel_size,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_candidates_produce_empty_matrix() {
        let ctx = FaCandidateContext::build(&[]);
        let no_reads: &[ReadWithPosInfo] = &[];
        let m = build_fa_matrix(&ctx, b"ACGT", 0, no_reads, &[], 89, 50);
        assert_eq!(m.matrix.len(), 0);
        assert_eq!(m.candidate_num, 0);
        assert!(m.alt_info_strings.is_empty());
    }
}
