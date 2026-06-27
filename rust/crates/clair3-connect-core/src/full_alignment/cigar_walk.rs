use super::candidates::FaCandidateContext;
use super::constants::normalize_bq;
use super::read_iter::RawRead;
use super::tables::{acgt_num_for_char, seqi_to_char};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Default)]
pub struct PosInfo {

    pub alt_base: i8,

    pub bq: i8,

    pub del_length: u32,

    pub ins_bases: Vec<u8>,

    pub signal_length: i32,
}

#[derive(Debug, Default, Clone)]
pub struct PosAltInfo {

    pub acgt_count: [u32; 4],

    pub depth: u32,

    pub ins_counter: BTreeMap<Vec<u8>, u32>,

    pub del_counter: BTreeMap<u32, u32>,
}

#[derive(Debug, Clone)]
pub struct ReadWithPosInfo {

    pub raw: RawRead,

    pub overlap_candidates_num: usize,

    pub flanking_start: usize,

    pub pos_info: Vec<PosInfo>,
}

#[derive(Debug, Clone)]
pub struct ProjectedRead {

    pub read_start: i64,

    pub read_end: i64,

    pub mq: i8,

    pub strand: i8,

    pub haplotype: u8,

    pub overlap_candidates_num: usize,

    pub flanking_start: usize,

    pub pos_info: Vec<PosInfo>,
}

pub fn process_read(
    read: RawRead,
    ctx: &FaCandidateContext,
    pos_alt_info: &mut [PosAltInfo],
) -> Option<ReadWithPosInfo> {
    let (flanking_start, overlap_candidates_num, pos_info) =
        cigar_walk_core(&read, ctx, pos_alt_info)?;
    Some(ReadWithPosInfo {
        raw: read,
        overlap_candidates_num,
        flanking_start,
        pos_info,
    })
}

pub fn process_read_projected(
    read: RawRead,
    ctx: &FaCandidateContext,
    pos_alt_info: &mut [PosAltInfo],
) -> Option<ProjectedRead> {
    let (flanking_start, overlap_candidates_num, pos_info) =
        cigar_walk_core(&read, ctx, pos_alt_info)?;
    Some(ProjectedRead {
        read_start: read.read_start,
        read_end: read.read_end,
        mq: read.mq,
        strand: read.strand,
        haplotype: read.haplotype,
        overlap_candidates_num,
        flanking_start,
        pos_info,
    })

}

fn cigar_walk_core(
    read: &RawRead,
    ctx: &FaCandidateContext,
    pos_alt_info: &mut [PosAltInfo],
) -> Option<(usize, usize, Vec<PosInfo>)> {
    let read_end = read.read_end as u64;
    let read_start = read.read_start as u64;

    let flanking_start = match ctx
        .flanking_positions
        .binary_search_by(|p| p.cmp(&read_start))
    {
        Ok(i) => i,
        Err(i) => i,
    };

    let mut overlap_candidates_num = 0;
    for i in flanking_start..ctx.flanking_positions.len() {
        let p = ctx.flanking_positions[i];
        if p >= read_start && p < read_end {
            overlap_candidates_num += 1;
        } else {
            break;
        }
    }
    if overlap_candidates_num == 0 {
        return None;
    }

    let mut pos_info: Vec<PosInfo> = vec![PosInfo::default(); overlap_candidates_num];

    let mut ref_pos: u64 = read_start;
    let mut query_pos: u32 = 0;

    for (op_char, length) in &read.cigar {
        let length = *length;
        let op = *op_char as u8 as char;
        match op {
            'M' | '=' | 'X' => {
                for p in ref_pos..ref_pos + length as u64 {
                    if let Some(flanking_index) = ctx.flanking_idx(p) {
                        if flanking_index >= flanking_start {
                            let offset = flanking_index - flanking_start;
                            if offset < overlap_candidates_num {
                                let seqi = bam_seqi(&read.seqi, query_pos);
                                pos_info[offset].alt_base = seqi as i8;
                                pos_info[offset].bq = if query_pos < read.qual.len() as u32 {
                                    normalize_bq(read.qual[query_pos as usize])
                                } else {
                                    0
                                };

                                if let Some(sigs) = &read.signal_lengths {
                                    if (query_pos as usize) < sigs.len() {
                                        pos_info[offset].signal_length =
                                            sigs[query_pos as usize];
                                    }
                                }
                                if let Some(c_idx) = ctx.candidate_idx(p) {
                                    let alt_char = seqi_to_char(seqi);
                                    let bucket = acgt_num_for_char(alt_char);
                                    if (0..4).contains(&bucket) {
                                        pos_alt_info[c_idx].acgt_count[bucket as usize] += 1;
                                    }
                                    pos_alt_info[c_idx].depth += 1;
                                }
                            }
                        }
                    }
                    query_pos += 1;
                }
                ref_pos += length as u64;
            }
            'D' => {

                if ref_pos >= 1 {
                    let anchor = ref_pos - 1;
                    if let Some(flanking_index) = ctx.flanking_idx(anchor) {
                        if flanking_index >= flanking_start {
                            let offset = flanking_index - flanking_start;
                            if offset < overlap_candidates_num {
                                pos_info[offset].del_length = length;
                                if let Some(c_idx) = ctx.candidate_idx(anchor) {
                                    *pos_alt_info[c_idx]
                                        .del_counter
                                        .entry(length)
                                        .or_insert(0) += 1;
                                }
                            }
                        }
                    }
                }

                for p in ref_pos..ref_pos + length as u64 {
                    if let Some(flanking_index) = ctx.flanking_idx(p) {
                        if flanking_index >= flanking_start {
                            let offset = flanking_index - flanking_start;
                            if offset < overlap_candidates_num {
                                pos_info[offset].alt_base = -1;
                                if let Some(c_idx) = ctx.candidate_idx(p) {
                                    pos_alt_info[c_idx].depth += 1;
                                }
                            }
                        }
                    }
                }
                ref_pos += length as u64;
            }
            'I' => {
                if ref_pos >= 1 {
                    let anchor = ref_pos - 1;
                    if let Some(flanking_index) = ctx.flanking_idx(anchor) {
                        if flanking_index >= flanking_start {
                            let offset = flanking_index - flanking_start;
                            if offset < overlap_candidates_num {
                                let mut ins_bases = Vec::with_capacity(length as usize);
                                let mut ins_signal_sum: i32 = 0;
                                for k in 0..length {
                                    let seqi = bam_seqi(&read.seqi, query_pos + k);
                                    ins_bases.push(seqi_to_char(seqi));

                                    if let Some(sigs) = &read.signal_lengths {
                                        let qp = (query_pos + k) as usize;
                                        if qp < sigs.len() {
                                            ins_signal_sum += sigs[qp];
                                        }
                                    }
                                }
                                pos_info[offset].ins_bases = ins_bases.clone();
                                if read.signal_lengths.is_some() {
                                    pos_info[offset].signal_length += ins_signal_sum;
                                }
                                if let Some(c_idx) = ctx.candidate_idx(anchor) {
                                    *pos_alt_info[c_idx]
                                        .ins_counter
                                        .entry(ins_bases)
                                        .or_insert(0) += 1;
                                }
                            }
                        }
                    }
                }
                query_pos += length;
            }
            'N' => {
                ref_pos += length as u64;
            }
            'S' => {
                query_pos += length;
            }
            'H' | 'P' => {

            }
            _ => {

            }
        }
    }

    Some((flanking_start, overlap_candidates_num, pos_info))
}

#[inline]
fn bam_seqi(seqi: &[u8], qpos: u32) -> u8 {
    let i = qpos as usize;
    let byte = seqi[i / 2];
    if i % 2 == 0 {
        byte >> 4
    } else {
        byte & 0x0F
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::full_alignment::read_iter::RawRead;

    fn rread(start: i64, end: i64, cigar: Vec<(u32, u32)>, seq: Vec<u8>) -> RawRead {
        RawRead {
            qname: b"r1".to_vec(),
            mq: 100,
            raw_mq: 60,
            read_start: start,
            read_end: end,
            cigar,
            seqi: seq,
            qual: vec![30; 200],
            l_qseq: 0,
            strand: 100,
            haplotype: 0,
            signal_lengths: None,
        }
    }

    fn pack_seq(bases: &[u8]) -> Vec<u8> {
        let table: [u8; 256] = {
            let mut t = [15u8; 256];
            t[b'A' as usize] = 1;
            t[b'C' as usize] = 2;
            t[b'G' as usize] = 4;
            t[b'T' as usize] = 8;
            t
        };
        let mut out = Vec::with_capacity(bases.len().div_ceil(2));
        for chunk in bases.chunks(2) {
            let hi = table[chunk[0] as usize];
            let lo = if chunk.len() == 2 { table[chunk[1] as usize] } else { 0 };
            out.push((hi << 4) | lo);
        }
        out
    }

    #[test]
    fn match_read_increments_acgt_at_candidate() {

        let ctx = FaCandidateContext::build(&[100]);
        let bases = b"AAAACATGCTA";
        let read = rread(
            95,
            106,
            vec![('M' as u32, 11)],
            pack_seq(bases),
        );
        let mut pai = vec![PosAltInfo::default(); 1];
        let r = process_read(read, &ctx, &mut pai).expect("overlap expected");

        assert_eq!(pai[0].acgt_count, [1, 0, 0, 0], "acgt: {:?}", pai[0]);
        assert_eq!(pai[0].depth, 1);

        assert!(r.overlap_candidates_num > 0);
    }

    #[test]
    fn deletion_anchored_at_prev_position() {

        let ctx = FaCandidateContext::build(&[99]);
        let bases: Vec<u8> = b"AAAACATGCC".iter().chain(b"GGGGGGG".iter()).copied().collect();
        let read = rread(
            90,
            110,
            vec![('M' as u32, 10), ('D' as u32, 3), ('M' as u32, 7)],
            pack_seq(&bases),
        );
        let mut pai = vec![PosAltInfo::default(); 1];
        let r = process_read(read, &ctx, &mut pai).expect("overlap");

        assert_eq!(pai[0].acgt_count[1], 1);
        assert_eq!(pai[0].del_counter.get(&3), Some(&1));

        assert_eq!(pai[0].depth, 1);
        assert!(r.overlap_candidates_num > 0);
    }

    #[test]
    fn insertion_anchored_at_prev_position() {

        let ctx = FaCandidateContext::build(&[99]);
        let bases: Vec<u8> = b"AAAACATGCC".iter().chain(b"TC".iter()).chain(b"GGGGG".iter()).copied().collect();
        let read = rread(
            90,
            105,
            vec![('M' as u32, 10), ('I' as u32, 2), ('M' as u32, 5)],
            pack_seq(&bases),
        );
        let mut pai = vec![PosAltInfo::default(); 1];
        let _ = process_read(read, &ctx, &mut pai).expect("overlap");

        assert_eq!(pai[0].acgt_count[1], 1);
        assert_eq!(pai[0].ins_counter.get(b"TC".as_slice()), Some(&1));
    }

    #[test]
    fn no_overlap_returns_none() {
        let ctx = FaCandidateContext::build(&[1000]);
        let bases = b"AAAAA";
        let read = rread(10, 15, vec![('M' as u32, 5)], pack_seq(bases));
        let mut pai = vec![PosAltInfo::default(); 1];
        assert!(process_read(read, &ctx, &mut pai).is_none());
    }

    #[test]
    fn soft_clip_advances_query_only() {

        let ctx = FaCandidateContext::build(&[100]);
        let bases = b"XXXAAAAA";
        let read = rread(
            100,
            105,
            vec![('S' as u32, 3), ('M' as u32, 5)],
            pack_seq(bases),
        );
        let mut pai = vec![PosAltInfo::default(); 1];
        let _ = process_read(read, &ctx, &mut pai).expect("overlap");

        assert_eq!(pai[0].acgt_count, [1, 0, 0, 0]);
    }

    #[test]
    fn ref_skip_advances_ref_only() {

        let ctx = FaCandidateContext::build(&[101]);
        let bases = b"AAAAATTTTT";
        let read = rread(
            95,
            108,
            vec![('M' as u32, 5), ('N' as u32, 3), ('M' as u32, 5)],
            pack_seq(bases),
        );
        let mut pai = vec![PosAltInfo::default(); 1];
        let _ = process_read(read, &ctx, &mut pai).expect("overlap");

        assert_eq!(pai[0].acgt_count, [0, 0, 0, 0]);
        assert_eq!(pai[0].depth, 0);
    }

    #[test]
    fn complex_cigar_nmnmd() {

        let ctx = FaCandidateContext::build(&[107]);
        let bases = b"AAACCCGG";
        let read = rread(
            100,
            111,
            vec![
                ('M' as u32, 3),
                ('N' as u32, 2),
                ('M' as u32, 3),
                ('D' as u32, 1),
                ('M' as u32, 2),
            ],
            pack_seq(bases),
        );
        let mut pai = vec![PosAltInfo::default(); 1];
        let _ = process_read(read, &ctx, &mut pai).expect("overlap");

        assert_eq!(pai[0].acgt_count[1], 1);
        assert_eq!(pai[0].del_counter.get(&1), Some(&1));
    }

    #[test]
    fn hard_clip_pad_are_noops() {

        let ctx = FaCandidateContext::build(&[100]);
        let bases = b"AAAAA";
        let read = rread(
            100,
            105,
            vec![('H' as u32, 2), ('M' as u32, 5), ('P' as u32, 3)],
            pack_seq(bases),
        );
        let mut pai = vec![PosAltInfo::default(); 1];
        let _ = process_read(read, &ctx, &mut pai).expect("overlap");

        assert_eq!(pai[0].acgt_count, [1, 0, 0, 0]);
    }
}
