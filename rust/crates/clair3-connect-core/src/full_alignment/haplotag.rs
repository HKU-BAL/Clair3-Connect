use super::constants::{HAP_1, HAP_2, HAP_UNPHASED, OVERHANG};
use super::levenshtein::levenshtein;
use super::phase_block_index::Variant;
use std::collections::HashMap;

mod op {
    pub const M: u32 = 0;
    pub const I: u32 = 1;
    pub const D: u32 = 2;
    pub const N: u32 = 3;
    pub const S: u32 = 4;
    pub const EQ: u32 = 7;
    pub const X: u32 = 8;
}

fn op_code_from_char(c: u32) -> u32 {
    match c as u8 as char {
        'M' => op::M,
        'I' => op::I,
        'D' => op::D,
        'N' => op::N,
        'S' => op::S,
        '=' => op::EQ,
        'X' => op::X,

        _ => op::M + 32,
    }
}

pub fn cigar_prefix_length(
    cigar: &[(u32, u32)],
    reference_bases: usize,
    left_cigar_index: usize,
    right_cigar_index: usize,
    consumed: usize,
    reverse: bool,
) -> (usize, usize) {
    let mut ref_pos: usize = 0;
    let mut query_pos: usize = 0;
    for i in left_cigar_index..right_cigar_index {
        let index = if reverse {
            left_cigar_index + right_cigar_index - i - 1
        } else {
            i
        };
        let (op_char, raw_len) = cigar[index];
        let cigar_op = op_code_from_char(op_char);
        let length = if i == left_cigar_index { consumed } else { raw_len as usize };
        if length == 0 {
            continue;
        }
        match cigar_op {
            op::M | op::EQ | op::X => {
                query_pos += length;
                ref_pos += length;
                if ref_pos >= reference_bases {
                    return (reference_bases, query_pos + reference_bases - ref_pos);
                }
            }
            op::D => {
                ref_pos += length;
                if ref_pos >= reference_bases {
                    return (reference_bases, query_pos);
                }
            }
            op::I => {
                query_pos += length;
            }
            op::N => {
                return (reference_bases, query_pos);
            }
            _ => {}
        }
    }
    (ref_pos, query_pos)
}

fn decode_query(seqi: &[u8], start: usize, end: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(end.saturating_sub(start));
    for q in start..end {
        let byte = seqi[q / 2];
        let nibble = if q % 2 == 0 { byte >> 4 } else { byte & 0x0F };
        out.push(super::tables::seqi_to_char(nibble));
    }
    out
}

#[allow(clippy::too_many_arguments)]
pub fn realign_read(
    variant: &Variant,
    cigar: &[(u32, u32)],
    seqi: &[u8],
    i: usize,
    consumed: usize,
    query_pos: usize,
    reference: &[u8],
    ref_start: i64,
) -> i32 {
    let n_cigar = cigar.len();
    let middle_length = cigar[i].1 as usize;
    let left_consumed = consumed.min(middle_length);
    let right_consumed = middle_length.saturating_sub(consumed);
    let left_cigar_size = i + 1;
    let right_cigar_size = i;

    let (left_ref_bases, left_query_bases) =
        cigar_prefix_length(cigar, OVERHANG, 0, left_cigar_size, left_consumed, true);
    let (right_ref_bases, right_query_bases) = cigar_prefix_length(
        cigar,
        OVERHANG + 1,
        right_cigar_size,
        n_cigar,
        right_consumed,
        false,
    );

    let qst = query_pos.saturating_sub(left_query_bases);
    let qen = query_pos + right_query_bases;
    if qen == qst {
        return 0;
    }
    let rst_signed = variant.position - left_ref_bases as i64 - ref_start;
    let ren_signed = variant.position + right_ref_bases as i64 - ref_start;
    if rst_signed < 0 || ren_signed < 0 || (ren_signed as usize) > reference.len() {
        return 0;
    }
    let rst = rst_signed as usize;
    let ren = ren_signed as usize;

    let query = decode_query(seqi, qst, qen);
    let r_ref = &reference[rst..ren];

    let alt_length = left_ref_bases + right_ref_bases;
    let alt_size = alt_length.max(left_ref_bases + 1);
    let mut alt: Vec<u8> = vec![0u8; alt_size];
    let copy_len = (ren - rst).min(alt_length);
    alt[..copy_len].copy_from_slice(&r_ref[..copy_len]);
    alt[left_ref_bases] = variant.alt_base;

    let distance_ref = levenshtein(&query, r_ref);
    let distance_alt = levenshtein(&query, &alt);

    if distance_ref < distance_alt {
        1
    } else if distance_ref > distance_alt {
        2
    } else {
        0
    }
}

fn update_haplotype_cost(allele: i32, phase_set: i32, genotype: i32, costs: &mut HashMap<i32, i32>) {
    if allele == 0 {
        return;
    }
    let entry = costs.entry(phase_set).or_insert(0);
    if allele == genotype {
        *entry += 1;
    } else {
        *entry -= 1;
    }
}

pub fn haplotag_read(
    variants: &[Variant],
    variant_current_pos: usize,
    read_start: i64,
    cigar: &[(u32, u32)],
    seqi: &[u8],
    reference: &[u8],
    ref_start: i64,
) -> (u8, i64) {
    let n = variants.len();
    let mut j = variant_current_pos;
    let mut ref_pos = read_start;
    let mut query_pos: usize = 0;
    let mut costs: HashMap<i32, i32> = HashMap::new();

    while j < n && variants[j].position < ref_pos {
        j += 1;
    }

    for i in 0..cigar.len() {
        let (op_char, length_u32) = cigar[i];
        let cigar_op = op_code_from_char(op_char);
        let length = length_u32 as i64;
        let mut v_position: i64 = if j < n { variants[j].position } else { i64::MAX };

        match cigar_op {
            op::M | op::EQ | op::X => {
                while j < n && v_position < ref_pos + length {
                    let v = &variants[j];
                    let local_offset = (v_position - ref_pos) as usize;
                    let allele = realign_read(
                        v,
                        cigar,
                        seqi,
                        i,
                        local_offset,
                        query_pos + local_offset,
                        reference,
                        ref_start,
                    );
                    update_haplotype_cost(allele, v.phase_set, v.genotype, &mut costs);
                    j += 1;
                    v_position = if j < n { variants[j].position } else { i64::MAX };
                }
                query_pos += length as usize;
                ref_pos += length;
            }
            op::I => {
                if j < n && v_position == ref_pos {
                    let v = &variants[j];
                    let allele = realign_read(v, cigar, seqi, i, 0, query_pos, reference, ref_start);
                    update_haplotype_cost(allele, v.phase_set, v.genotype, &mut costs);
                    j += 1;
                }
                query_pos += length as usize;
            }
            op::D => {
                while j < n && v_position < ref_pos + length {
                    let v = &variants[j];
                    let local_offset = (v_position - ref_pos) as usize;
                    let allele = realign_read(
                        v,
                        cigar,
                        seqi,
                        i,
                        local_offset,
                        query_pos,
                        reference,
                        ref_start,
                    );
                    update_haplotype_cost(allele, v.phase_set, v.genotype, &mut costs);
                    j += 1;
                    v_position = if j < n { variants[j].position } else { i64::MAX };
                }
                ref_pos += length;
            }
            op::N => {
                while j < n && v_position < ref_pos + length {
                    j += 1;
                    v_position = if j < n { variants[j].position } else { i64::MAX };
                }
                ref_pos += length;
            }
            op::S => {
                query_pos += length as usize;
            }
            _ => {}
        }
    }

    let (mut max_v, mut min_v) = (0_i32, 0_i32);
    let mut counter_size = 0;
    for v in costs.values() {
        max_v = max_v.max(*v);
        min_v = min_v.min(*v);
        counter_size += 1;
    }
    let hap = if counter_size == 0 || (max_v == 0 && min_v == 0) {
        HAP_UNPHASED
    } else if max_v > min_v.abs() {
        HAP_1
    } else {
        HAP_2
    };
    (hap, ref_pos)
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn cigar_prefix_simple_forward() {

        let cigar = vec![('M' as u32, 5u32)];
        let (r, q) = cigar_prefix_length(&cigar, 10, 0, 1, 5, false);
        assert_eq!((r, q), (5, 5));
    }

    #[test]
    fn cigar_prefix_simple_reverse() {

        let cigar = vec![('M' as u32, 5u32), ('I' as u32, 3), ('M' as u32, 5)];
        let (r, q) = cigar_prefix_length(&cigar, 4, 0, 3, 5, true);
        assert_eq!(r, 4);

        assert_eq!(q, 4);
    }

    #[test]
    fn haplotag_no_variants_returns_unphased() {
        let cigar = vec![('M' as u32, 20)];
        let seqi = pack_seq(b"AAAAAAAAAAAAAAAAAAAA");
        let r_ref = b"AAAAAAAAAAAAAAAAAAAA".as_slice();
        let (hap, _) = haplotag_read(&[], 0, 100, &cigar, &seqi, r_ref, 100);
        assert_eq!(hap, HAP_UNPHASED);
    }

    #[test]
    fn haplotag_one_snp_supports_alt_returns_hap1() {

        let variants = vec![Variant {
            position: 100,
            ref_base: b'A',
            alt_base: b'G',
            genotype: 1,
            phase_set: 200,
        }];
        let cigar = vec![('M' as u32, 21u32)];

        let query_bases = b"GAAAAAAAAAAAAAAAAAAAA";
        let seqi = pack_seq(query_bases);

        let r_ref: Vec<u8> = vec![b'A'; 35];
        let (hap, _) = haplotag_read(&variants, 0, 100, &cigar, &seqi, &r_ref, 90);

        assert_eq!(hap, HAP_2);
    }

    #[test]
    fn haplotag_one_snp_supports_ref_returns_hap1_or_hap2() {

        let variants = vec![Variant {
            position: 100,
            ref_base: b'A',
            alt_base: b'G',
            genotype: 1,
            phase_set: 200,
        }];
        let cigar = vec![('M' as u32, 21u32)];
        let query_bases = b"AAAAAAAAAAAAAAAAAAAAA";
        let seqi = pack_seq(query_bases);
        let r_ref: Vec<u8> = vec![b'A'; 35];
        let (hap, _) = haplotag_read(&variants, 0, 100, &cigar, &seqi, &r_ref, 90);
        assert_eq!(hap, HAP_1);
    }
}
