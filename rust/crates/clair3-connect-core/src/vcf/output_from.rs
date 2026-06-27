use super::find_alt_base::{find_alt_base, AltInfoEntries};
use super::indel_bases::{
    deletion_bases_using_alt_info_from, insertion_bases_using_alt_info_from, DeletionOptions,
    InsertionLookupResult, InsertionOptions,
};
use super::iupac::iupac_to_acgt;
use super::labels::variant_length;
use super::probabilities::{possible_outcomes_fa, possible_outcomes_pileup, FullOutcome, Outcome};
use super::snp::{hetero_snp_bases_from, homo_snp_bases_from};

#[derive(Debug, Clone, PartialEq)]
pub struct OutputDecision {

    pub is_reference: bool,

    pub is_homo_snp: bool,

    pub is_hetero_snp: bool,

    pub is_homo_insertion: bool,

    pub is_hetero_acgt_ins: bool,

    pub is_hetero_insins: bool,

    pub is_homo_deletion: bool,

    pub is_hetero_acgt_del: bool,

    pub is_hetero_deldel: bool,

    pub is_insertion_and_deletion: bool,

    pub reference_base: String,

    pub alternate_base: String,

    pub maximum_probability: f32,
}

impl OutputDecision {

    fn ref_call(ref_base: char, p: f32) -> Self {
        Self {
            is_reference: true,
            is_homo_snp: false,
            is_hetero_snp: false,
            is_homo_insertion: false,
            is_hetero_acgt_ins: false,
            is_hetero_insins: false,
            is_homo_deletion: false,
            is_hetero_acgt_del: false,
            is_hetero_deldel: false,
            is_insertion_and_deletion: false,
            reference_base: ref_base.to_string(),
            alternate_base: ref_base.to_string(),
            maximum_probability: p,
        }
    }
}

pub fn output_from(
    reference_sequence: &str,
    tensor_position_center: usize,
    gt21_probabilities: &[f32],
    genotype_probabilities: &[f32],
    entries: &AltInfoEntries,
    maximum_variant_length_that_need_infer: usize,
) -> OutputDecision {
    output_from_inner(
        reference_sequence,
        tensor_position_center,
        gt21_probabilities,
        genotype_probabilities,
        None,
        entries,
        maximum_variant_length_that_need_infer,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn output_from_fa(
    reference_sequence: &str,
    tensor_position_center: usize,
    gt21_probabilities: &[f32],
    genotype_probabilities: &[f32],
    vlen1: &[f32],
    vlen2: &[f32],
    entries: &AltInfoEntries,
    maximum_variant_length_that_need_infer: usize,
) -> OutputDecision {
    output_from_inner(
        reference_sequence,
        tensor_position_center,
        gt21_probabilities,
        genotype_probabilities,
        Some((vlen1, vlen2)),
        entries,
        maximum_variant_length_that_need_infer,
    )
}

#[allow(clippy::too_many_arguments)]
fn output_from_inner(
    reference_sequence: &str,
    tensor_position_center: usize,
    gt21_probabilities: &[f32],
    genotype_probabilities: &[f32],
    vlen_heads: Option<(&[f32], &[f32])>,
    entries: &AltInfoEntries,
    maximum_variant_length_that_need_infer: usize,
) -> OutputDecision {
    let ref_byte = reference_sequence.as_bytes()[tensor_position_center];

    debug_assert!(
        matches!(ref_byte, b'A' | b'C' | b'G' | b'T' | b'U' | b'R' | b'Y' | b'S' | b'W'
            | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' | b'N'),
        "output_from: reference byte {:?} is not in IUPAC; upstream caller \
         must pre-filter (Python BASIC_BASES); see CallVariants.py:1053",
        ref_byte as char
    );
    let ref_base_acgt_byte = iupac_to_acgt(ref_byte);

    let add_indel_length = vlen_heads.is_some();
    let outcome = if let Some((vlen1, vlen2)) = vlen_heads {
        possible_outcomes_fa(
            gt21_probabilities,
            genotype_probabilities,
            vlen1,
            vlen2,
            ref_base_acgt_byte,
        )
    } else {
        possible_outcomes_pileup(
            gt21_probabilities,
            genotype_probabilities,
            ref_base_acgt_byte,
        )
    };

    let mut full = match outcome {
        Outcome::HomoRefShortcut(p) => {
            return OutputDecision::ref_call(ref_base_acgt_byte as char, p);
        }
        Outcome::Full(f) => *f,
    };

    loop {
        let max_probability = compute_max(&full);
        let is_reference = max_probability == full.homo_ref;
        if is_reference {
            return OutputDecision::ref_call(ref_base_acgt_byte as char, max_probability);
        }

        let is_homo_snp = full.homo_snp.contains(&max_probability);
        let is_hetero_snp = full.hetero_snp.contains(&max_probability);
        let is_homo_insertion = full.homo_ins_probabilities.contains(&max_probability);
        let is_hetero_acgt_ins =
            full.hetero_acgt_ins_probabilities.contains(&max_probability);
        let is_hetero_insins =
            full.hetero_insins_probabilities.contains(&max_probability);
        let is_homo_deletion = full.homo_del_probabilities.contains(&max_probability);
        let is_hetero_acgt_del =
            full.hetero_acgt_del_probabilities.contains(&max_probability);
        let is_hetero_deldel =
            full.hetero_deldel_probabilities.contains(&max_probability);
        let is_insertion_and_deletion =
            full.hetero_insdel_probabilities.contains(&max_probability);

        let raw_ref_char = ref_byte as char;

        if is_homo_snp {
            let idx = first_index(&full.homo_snp, max_probability);
            let (base1, base2) = homo_snp_bases_from(&full.homo_snp);
            let alt = if base1 != raw_ref_char { base1 } else { base2 };
            let (_sorted, rewritten) = find_alt_base(entries, Some(alt));
            match rewritten {
                Some(rb) if rb != raw_ref_char => {
                    return decide(
                        is_homo_snp,
                        is_hetero_snp,
                        is_homo_insertion,
                        is_hetero_acgt_ins,
                        is_hetero_insins,
                        is_homo_deletion,
                        is_hetero_acgt_del,
                        is_hetero_deldel,
                        is_insertion_and_deletion,
                        raw_ref_char.to_string(),
                        rb.to_string(),
                        max_probability,
                    );
                }
                _ => {
                    full.homo_snp[idx] = 0.0;
                    continue;
                }
            }
        }

        if is_hetero_snp {
            let idx = first_index(&full.hetero_snp, max_probability);
            let (base1, base2) = hetero_snp_bases_from(&full.hetero_snp);
            let is_multi = base1 != raw_ref_char && base2 != raw_ref_char;
            if is_multi {
                let (sorted_alt_bases, _) = find_alt_base(entries, None);
                if sorted_alt_bases.len() < 2 {
                    full.hetero_snp[idx] = 0.0;
                    continue;
                }
                let alt = format!("{},{}", sorted_alt_bases[0], sorted_alt_bases[1]);
                return decide(
                    is_homo_snp,
                    is_hetero_snp,
                    is_homo_insertion,
                    is_hetero_acgt_ins,
                    is_hetero_insins,
                    is_homo_deletion,
                    is_hetero_acgt_del,
                    is_hetero_deldel,
                    is_insertion_and_deletion,
                    raw_ref_char.to_string(),
                    alt,
                    max_probability,
                );
            }
            let alt = if base1 != raw_ref_char { base1 } else { base2 };
            let (_sorted, rewritten) = find_alt_base(entries, Some(alt));
            match rewritten {
                Some(rb) if rb != raw_ref_char => {
                    return decide(
                        is_homo_snp,
                        is_hetero_snp,
                        is_homo_insertion,
                        is_hetero_acgt_ins,
                        is_hetero_insins,
                        is_homo_deletion,
                        is_hetero_acgt_del,
                        is_hetero_deldel,
                        is_insertion_and_deletion,
                        raw_ref_char.to_string(),
                        rb.to_string(),
                        max_probability,
                    );
                }
                _ => {
                    full.hetero_snp[idx] = 0.0;
                    continue;
                }
            }
        }

        if is_homo_insertion {
            let idx = first_index_slice(&full.homo_ins_probabilities, max_probability);

            let propose_len = if add_indel_length && idx < full.homo_ins_lengths.len() {
                let l = full.homo_ins_lengths[idx];
                if l > 0 && l < variant_length::MAX {
                    Some(l as usize)
                } else {
                    None
                }
            } else {
                None
            };
            let ins = insertion_bases_using_alt_info_from(
                entries,
                InsertionOptions {
                    propose_insertion_length: propose_len,
                    maximum_insertion_length: maximum_variant_length_that_need_infer,
                    ..Default::default()
                },
            );
            if ins.is_empty() {
                full.homo_ins_probabilities[idx] = 0.0;
                continue;
            }
            return decide(
                is_homo_snp,
                is_hetero_snp,
                is_homo_insertion,
                is_hetero_acgt_ins,
                is_hetero_insins,
                is_homo_deletion,
                is_hetero_acgt_del,
                is_hetero_deldel,
                is_insertion_and_deletion,
                raw_ref_char.to_string(),
                ins.best().to_string(),
                max_probability,
            );
        }

        if is_hetero_acgt_ins {
            let idx = first_index_slice(&full.hetero_acgt_ins_probabilities, max_probability);

            let hetero_ins_base = if add_indel_length && idx < full.hetero_acgt_ins_bases.len() {
                full.hetero_acgt_ins_bases[idx]
            } else {
                ACGT[idx % 4]
            };
            let propose_len = if add_indel_length && idx < full.hetero_acgt_ins_lengths.len() {
                let l = full.hetero_acgt_ins_lengths[idx];
                if l > 0 && l < variant_length::MAX {
                    Some(l as usize)
                } else {
                    None
                }
            } else {
                None
            };
            let ins = insertion_bases_using_alt_info_from(
                entries,
                InsertionOptions {
                    propose_insertion_length: propose_len,
                    maximum_insertion_length: maximum_variant_length_that_need_infer,
                    ..Default::default()
                },
            );
            if ins.is_empty() {
                full.hetero_acgt_ins_probabilities[idx] = 0.0;
                continue;
            }
            let mut alt = ins.best().to_string();

            let is_snp_ins_multi = hetero_ins_base != raw_ref_char;
            if is_snp_ins_multi {
                let (sorted_alt_bases, _) = find_alt_base(entries, None);
                if sorted_alt_bases.is_empty() {
                    full.hetero_acgt_ins_probabilities[idx] = 0.0;
                    continue;
                }
                alt = format!("{},{}", sorted_alt_bases[0], alt);
            }
            return decide(
                is_homo_snp,
                is_hetero_snp,
                is_homo_insertion,
                is_hetero_acgt_ins,
                is_hetero_insins,
                is_homo_deletion,
                is_hetero_acgt_del,
                is_hetero_deldel,
                is_insertion_and_deletion,
                raw_ref_char.to_string(),
                alt,
                max_probability,
            );
        }

        if is_hetero_insins {
            let idx = first_index_slice(&full.hetero_insins_probabilities, max_probability);

            let mut ins_list: Vec<String> = Vec::new();
            let mut fa_resolved = false;
            if add_indel_length && idx < full.hetero_insins_length_tuples.len() {
                let (l1, l2) = full.hetero_insins_length_tuples[idx];
                let p1 = if l1 > 0 && l1 < variant_length::MAX {
                    Some(l1 as usize)
                } else {
                    None
                };
                let p2 = if l2 > 0 && l2 < variant_length::MAX {
                    Some(l2 as usize)
                } else {
                    None
                };
                let ins1 = insertion_bases_using_alt_info_from(
                    entries,
                    InsertionOptions {
                        propose_insertion_length: p1,
                        maximum_insertion_length: maximum_variant_length_that_need_infer,
                        ..Default::default()
                    },
                );
                if !ins1.is_empty() {
                    let bases1 = ins1.best().to_string();
                    let ins2 = insertion_bases_using_alt_info_from(
                        entries,
                        InsertionOptions {
                            propose_insertion_length: p2,
                            insertion_bases_to_ignore: &bases1,
                            maximum_insertion_length: maximum_variant_length_that_need_infer,
                            ..Default::default()
                        },
                    );
                    if !ins2.is_empty() {
                        ins_list = vec![bases1, ins2.best().to_string()];
                        fa_resolved = true;
                    }
                }
            }
            if !fa_resolved {
                let multi = insertion_bases_using_alt_info_from(
                    entries,
                    InsertionOptions {
                        return_multi: true,
                        maximum_insertion_length: maximum_variant_length_that_need_infer,
                        ..Default::default()
                    },
                );
                ins_list = match multi {
                    InsertionLookupResult::Multi(v) => v,
                    _ => Vec::new(),
                };
            }
            if ins_list.len() < 2 {
                full.hetero_insins_probabilities[idx] = 0.0;
                continue;
            }

            let insertion_bases = &ins_list[0];
            let another = &ins_list[1];
            if insertion_bases == another {
                full.hetero_insins_probabilities[idx] = 0.0;
                continue;
            }

            let alt = format!("{},{}", another, insertion_bases);
            return decide(
                is_homo_snp,
                is_hetero_snp,
                is_homo_insertion,
                is_hetero_acgt_ins,
                is_hetero_insins,
                is_homo_deletion,
                is_hetero_acgt_del,
                is_hetero_deldel,
                is_insertion_and_deletion,
                raw_ref_char.to_string(),
                alt,
                max_probability,
            );
        }

        if is_homo_deletion {
            let idx = first_index_slice(&full.homo_del_probabilities, max_probability);
            let propose_len = if add_indel_length && idx < full.homo_del_lengths.len() {
                let l = full.homo_del_lengths[idx];
                if l > 0 && l < variant_length::MAX {
                    Some(l as usize)
                } else {
                    None
                }
            } else {
                None
            };
            let del = deletion_bases_using_alt_info_from(
                entries,
                DeletionOptions {
                    propose_deletion_length: propose_len,
                    maximum_deletion_length: maximum_variant_length_that_need_infer,
                    ..Default::default()
                },
            );
            let del_bases = del.best();
            if del_bases.is_empty() {
                full.homo_del_probabilities[idx] = 0.0;
                continue;
            }
            let mut reference_base = String::with_capacity(1 + del_bases.len());
            reference_base.push(raw_ref_char);
            reference_base.push_str(del_bases);
            let alternate_base = reference_base[..1].to_string();
            return decide(
                is_homo_snp,
                is_hetero_snp,
                is_homo_insertion,
                is_hetero_acgt_ins,
                is_hetero_insins,
                is_homo_deletion,
                is_hetero_acgt_del,
                is_hetero_deldel,
                is_insertion_and_deletion,
                reference_base,
                alternate_base,
                max_probability,
            );
        }

        if is_hetero_acgt_del {
            let idx = first_index_slice(&full.hetero_acgt_del_probabilities, max_probability);
            let hetero_del_base = if add_indel_length && idx < full.hetero_acgt_del_bases.len() {
                full.hetero_acgt_del_bases[idx]
            } else {
                ACGT[idx % 4]
            };
            let propose_len = if add_indel_length && idx < full.hetero_acgt_del_lengths.len() {
                let l = full.hetero_acgt_del_lengths[idx];
                if l > 0 && l < variant_length::MAX {
                    Some(l as usize)
                } else {
                    None
                }
            } else {
                None
            };
            let del = deletion_bases_using_alt_info_from(
                entries,
                DeletionOptions {
                    propose_deletion_length: propose_len,
                    maximum_deletion_length: maximum_variant_length_that_need_infer,
                    ..Default::default()
                },
            );
            let del_bases = del.best();
            if del_bases.is_empty() {
                full.hetero_acgt_del_probabilities[idx] = 0.0;
                continue;
            }
            let mut reference_base = String::with_capacity(1 + del_bases.len());
            reference_base.push(raw_ref_char);
            reference_base.push_str(del_bases);
            let mut alternate_base = reference_base[..1].to_string();

            let is_snp_del_multi = hetero_del_base != raw_ref_char;
            if is_snp_del_multi {

                let alt2 = format!("{}{}", hetero_del_base, &reference_base[1..]);
                alternate_base = format!("{},{}", alternate_base, alt2);
            }
            return decide(
                is_homo_snp,
                is_hetero_snp,
                is_homo_insertion,
                is_hetero_acgt_ins,
                is_hetero_insins,
                is_homo_deletion,
                is_hetero_acgt_del,
                is_hetero_deldel,
                is_insertion_and_deletion,
                reference_base,
                alternate_base,
                max_probability,
            );
        }

        if is_hetero_deldel {
            let idx = first_index_slice(&full.hetero_deldel_probabilities, max_probability);

            let mut del_list: Vec<String> = Vec::new();
            let mut fa_resolved = false;
            if add_indel_length && idx < full.hetero_deldel_length_tuples.len() {
                let (a, b) = full.hetero_deldel_length_tuples[idx];
                let (l1, l2) = if a > b { (a, b) } else { (b, a) };
                let p1 = if l1 > 0 && l1 < variant_length::MAX {
                    Some(l1 as usize)
                } else {
                    None
                };
                let p2 = if l2 > 0 && l2 < variant_length::MAX {
                    Some(l2 as usize)
                } else {
                    None
                };
                let del1 = deletion_bases_using_alt_info_from(
                    entries,
                    DeletionOptions {
                        propose_deletion_length: p1,
                        maximum_deletion_length: maximum_variant_length_that_need_infer,
                        ..Default::default()
                    },
                );
                let del1_bases = del1.best().to_string();
                if !del1_bases.is_empty() {
                    let del2 = deletion_bases_using_alt_info_from(
                        entries,
                        DeletionOptions {
                            propose_deletion_length: p2,
                            deletion_bases_to_ignore: &del1_bases,
                            maximum_deletion_length: maximum_variant_length_that_need_infer,
                            ..Default::default()
                        },
                    );
                    let del2_bases = del2.best().to_string();
                    if !del2_bases.is_empty() {

                        del_list = if del1_bases.len() > del2_bases.len() {
                            vec![del1_bases, del2_bases]
                        } else {
                            vec![del2_bases, del1_bases]
                        };
                        fa_resolved = true;
                    }
                }
            }
            if !fa_resolved {
                let multi = deletion_bases_using_alt_info_from(
                    entries,
                    DeletionOptions {
                        return_multi: true,
                        maximum_deletion_length: maximum_variant_length_that_need_infer,
                        ..Default::default()
                    },
                );
                del_list = match multi {
                    InsertionLookupResult::Multi(v) => v,
                    _ => Vec::new(),
                };
            }
            if del_list.len() < 2 {
                full.hetero_deldel_probabilities[idx] = 0.0;
                continue;
            }

            let deletion_bases = &del_list[0];
            let deletion_bases1 = &del_list[1];

            let mut reference_base = String::with_capacity(1 + deletion_bases.len());
            reference_base.push(raw_ref_char);
            reference_base.push_str(deletion_bases);

            let alt1 = reference_base[..1].to_string();

            let cut = deletion_bases1.len() + 1;
            let alt2 = if cut < reference_base.len() {
                format!("{}{}", &reference_base[..1], &reference_base[cut..])
            } else {
                reference_base[..1].to_string()
            };

            if alt1 == alt2 || reference_base == alt1 || reference_base == alt2 {
                full.hetero_deldel_probabilities[idx] = 0.0;
                continue;
            }
            let alternate_base = format!("{},{}", alt1, alt2);
            return decide(
                is_homo_snp,
                is_hetero_snp,
                is_homo_insertion,
                is_hetero_acgt_ins,
                is_hetero_insins,
                is_homo_deletion,
                is_hetero_acgt_del,
                is_hetero_deldel,
                is_insertion_and_deletion,
                reference_base,
                alternate_base,
                max_probability,
            );
        }

        if is_insertion_and_deletion {
            let idx = first_index_slice(&full.hetero_insdel_probabilities, max_probability);

            let (del_len_proposal, ins_len_proposal) = if add_indel_length
                && idx < full.hetero_insdel_length_tuples.len()
            {
                let (l1, l2) = full.hetero_insdel_length_tuples[idx];
                let dp = if l1 > 0 && l1 < variant_length::MAX {
                    Some(l1 as usize)
                } else {
                    None
                };
                let ip = if l2 > 0 && l2 < variant_length::MAX {
                    Some(l2 as usize)
                } else {
                    None
                };
                (dp, ip)
            } else {
                (None, None)
            };
            let ins = insertion_bases_using_alt_info_from(
                entries,
                InsertionOptions {
                    propose_insertion_length: ins_len_proposal,
                    maximum_insertion_length: maximum_variant_length_that_need_infer,
                    ..Default::default()
                },
            );
            let del = deletion_bases_using_alt_info_from(
                entries,
                DeletionOptions {
                    propose_deletion_length: del_len_proposal,
                    maximum_deletion_length: maximum_variant_length_that_need_infer,
                    ..Default::default()
                },
            );
            let ins_bases = ins.best();
            let del_bases = del.best();
            if ins_bases.is_empty() || del_bases.is_empty() {
                full.hetero_insdel_probabilities[idx] = 0.0;
                continue;
            }
            let mut reference_base = String::with_capacity(1 + del_bases.len());
            reference_base.push(raw_ref_char);
            reference_base.push_str(del_bases);

            let alt = format!(
                "{},{}{}",
                &reference_base[..1],
                ins_bases,
                &reference_base[1..]
            );
            return decide(
                is_homo_snp,
                is_hetero_snp,
                is_homo_insertion,
                is_hetero_acgt_ins,
                is_hetero_insins,
                is_homo_deletion,
                is_hetero_acgt_del,
                is_hetero_deldel,
                is_insertion_and_deletion,
                reference_base,
                alt,
                max_probability,
            );
        }

        return OutputDecision::ref_call(ref_base_acgt_byte as char, max_probability);
    }
}

const ACGT: [char; 4] = ['A', 'C', 'G', 'T'];

fn compute_max(f: &FullOutcome) -> f32 {

    let mut m = f.homo_ref;
    fn upd(m: &mut f32, v: f32) {
        if v > *m {
            *m = v;
        }
    }
    for v in &f.homo_snp {
        upd(&mut m, *v);
    }
    for v in &f.hetero_snp {
        upd(&mut m, *v);
    }
    for v in &f.homo_ins_probabilities {
        upd(&mut m, *v);
    }
    for v in &f.homo_del_probabilities {
        upd(&mut m, *v);
    }
    for v in &f.hetero_acgt_ins_probabilities {
        upd(&mut m, *v);
    }
    for v in &f.hetero_insins_probabilities {
        upd(&mut m, *v);
    }
    for v in &f.hetero_acgt_del_probabilities {
        upd(&mut m, *v);
    }
    for v in &f.hetero_deldel_probabilities {
        upd(&mut m, *v);
    }
    for v in &f.hetero_insdel_probabilities {
        upd(&mut m, *v);
    }
    m
}

fn first_index<const N: usize>(arr: &[f32; N], target: f32) -> usize {
    arr.iter().position(|x| *x == target).unwrap_or(0)
}

fn first_index_slice(arr: &[f32], target: f32) -> usize {
    arr.iter().position(|x| *x == target).unwrap_or(0)
}

#[allow(clippy::too_many_arguments)]
fn decide(
    is_homo_snp: bool,
    is_hetero_snp: bool,
    is_homo_insertion: bool,
    is_hetero_acgt_ins: bool,
    is_hetero_insins: bool,
    is_homo_deletion: bool,
    is_hetero_acgt_del: bool,
    is_hetero_deldel: bool,
    is_insertion_and_deletion: bool,
    reference_base: String,
    alternate_base: String,
    p: f32,
) -> OutputDecision {
    OutputDecision {
        is_reference: false,
        is_homo_snp,
        is_hetero_snp,
        is_homo_insertion,
        is_hetero_acgt_ins,
        is_hetero_insins,
        is_homo_deletion,
        is_hetero_acgt_del,
        is_hetero_deldel,
        is_insertion_and_deletion,
        reference_base,
        alternate_base,
        maximum_probability: p,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entries(pairs: &[(&str, u32)]) -> AltInfoEntries {
        pairs.iter().map(|(k, v)| (k.to_string(), *v)).collect()
    }

    fn one_hot_gt21(idx: usize, weight: f32) -> Vec<f32> {
        let mut v = vec![(1.0 - weight) / 20.0; 21];
        v[idx] = weight;
        v
    }

    fn one_hot_geno(idx: usize, weight: f32) -> Vec<f32> {
        let mut v = vec![(1.0 - weight) / 2.0; 3];
        v[idx] = weight;
        v
    }

    #[test]
    fn refcall_via_shortcut() {

        let gt21 = one_hot_gt21(0, 0.95);
        let geno = one_hot_geno(0, 0.9);
        let d = output_from("A", 0, &gt21, &geno, &entries(&[("RA", 100)]), 50);
        assert!(d.is_reference);
        assert_eq!(d.reference_base, "A");
        assert_eq!(d.alternate_base, "A");
    }

    #[test]
    fn homo_snp_resolves_via_alt_info() {

        use super::super::labels::Gt21;
        let gt21 = one_hot_gt21(Gt21::CC as usize, 0.9);
        let geno = one_hot_geno(1, 0.95);
        let d = output_from(
            "A",
            0,
            &gt21,
            &geno,
            &entries(&[("XC", 30), ("RA", 1)]),
            50,
        );
        assert!(d.is_homo_snp, "got {d:?}");
        assert_eq!(d.reference_base, "A");
        assert_eq!(d.alternate_base, "C");
    }

    #[test]
    fn hetero_snp_resolves_via_alt_info() {

        use super::super::labels::Gt21;
        let gt21 = one_hot_gt21(Gt21::AC as usize, 0.9);
        let geno = one_hot_geno(2, 0.95);
        let d = output_from(
            "A",
            0,
            &gt21,
            &geno,
            &entries(&[("XC", 20), ("RA", 10)]),
            50,
        );
        assert!(d.is_hetero_snp);
        assert_eq!(d.reference_base, "A");
        assert_eq!(d.alternate_base, "C");
    }

    #[test]
    fn hetero_snp_multi_allelic() {

        use super::super::labels::Gt21;
        let gt21 = one_hot_gt21(Gt21::CT as usize, 0.9);
        let geno = one_hot_geno(2, 0.95);
        let d = output_from(
            "A",
            0,
            &gt21,
            &geno,
            &entries(&[("XC", 20), ("XT", 15)]),
            50,
        );
        assert!(d.is_hetero_snp);
        assert_eq!(d.reference_base, "A");

        assert_eq!(d.alternate_base, "C,T");
    }

    #[test]
    fn homo_insertion_resolves() {
        use super::super::labels::Gt21;
        let gt21 = one_hot_gt21(Gt21::InsIns as usize, 0.85);
        let geno = one_hot_geno(1, 0.9);

        let d = output_from(
            "A",
            0,
            &gt21,
            &geno,
            &entries(&[("IAAT", 12), ("RA", 3)]),
            50,
        );
        assert!(d.is_homo_insertion);
        assert_eq!(d.reference_base, "A");
        assert_eq!(d.alternate_base, "AAT");
    }

    #[test]
    fn homo_deletion_resolves() {
        use super::super::labels::Gt21;
        let gt21 = one_hot_gt21(Gt21::DelDel as usize, 0.9);
        let geno = one_hot_geno(1, 0.95);

        let d = output_from(
            "A",
            0,
            &gt21,
            &geno,
            &entries(&[("DT", 15), ("RA", 2)]),
            50,
        );
        assert!(d.is_homo_deletion);
        assert_eq!(d.reference_base, "AT");
        assert_eq!(d.alternate_base, "A");
    }

    #[test]
    fn unresolvable_falls_back_to_ref_call() {

        use super::super::labels::Gt21;
        let gt21 = one_hot_gt21(Gt21::AT as usize, 0.7);
        let geno = one_hot_geno(2, 0.6);
        let d = output_from("A", 0, &gt21, &geno, &entries(&[("RA", 100)]), 50);

        assert!(d.is_reference || d.is_homo_snp || d.is_hetero_snp,
            "expected resolution: {d:?}");
    }
}
