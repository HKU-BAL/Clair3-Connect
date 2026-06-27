use super::find_alt_base::{decode_alt_info, AltInfoEntries, AltInfoError};
use super::iupac::convert_iupac_to_n;
use super::output_from::{output_from, output_from_fa, OutputDecision};
use super::quality::{filtration_value_from, quality_score_from};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum OutputWithError {

    #[error(transparent)]
    AltInfo(#[from] AltInfoError),

    #[error("chr_pos_seq malformed (need at least chr:pos:ref): {0:?}")]
    BadChrPosSeq(String),

    #[error("chr_pos_seq position not a u64: {0:?}")]
    BadPosition(String),
}

#[derive(Debug, Clone)]
pub struct OutputWithInputs<'a> {

    pub chr_pos_seq: &'a str,

    pub alt_info: &'a str,

    pub gt21_probabilities: &'a [f32],

    pub genotype_probabilities: &'a [f32],

    pub vlen_heads: Option<(&'a [f32], &'a [f32])>,
}

#[derive(Debug, Clone, Copy)]
pub struct OutputWithConfig {

    pub is_show_reference: bool,

    pub pileup: bool,

    pub add_indel_length: bool,

    pub quality_score_for_pass: Option<f64>,

    pub keep_iupac_bases: bool,

    pub is_haploid_precise_mode_enabled: bool,

    pub is_haploid_sensitive_mode_enabled: bool,

    pub enable_long_indel: bool,

    pub maximum_variant_length_that_need_infer: usize,
}

impl Default for OutputWithConfig {
    fn default() -> Self {
        Self {
            is_show_reference: false,
            pileup: true,
            add_indel_length: false,
            quality_score_for_pass: None,
            keep_iupac_bases: false,
            is_haploid_precise_mode_enabled: false,
            is_haploid_sensitive_mode_enabled: false,
            enable_long_indel: false,
            maximum_variant_length_that_need_infer: 50,
        }
    }
}

pub fn output_with(
    inputs: OutputWithInputs<'_>,
    config: &OutputWithConfig,
) -> Result<Option<String>, OutputWithError> {

    let info_list: Vec<&str> = inputs.chr_pos_seq.trim_end().split(':').collect();
    if info_list.len() < 3 {
        return Err(OutputWithError::BadChrPosSeq(
            inputs.chr_pos_seq.to_string(),
        ));
    }
    let (chromosome, position_str, reference_sequence) = if info_list.len() == 3 {
        (info_list[0].to_string(), info_list[1], info_list[2])
    } else {
        let n = info_list.len();

        (info_list[..n - 2].join(":"), info_list[n - 2], info_list[n - 1])
    };
    let position: u64 = position_str
        .parse()
        .map_err(|_| OutputWithError::BadPosition(position_str.to_string()))?;

    let (read_depth, entries) = decode_alt_info(inputs.alt_info.trim_end())?;

    let tensor_position_center = if reference_sequence.len() > 1 {

        16
    } else {
        0
    };
    let info_flag = if config.pileup { "P" } else { "F" };

    Ok(finish(
        &chromosome,
        position,
        reference_sequence,
        tensor_position_center,
        info_flag,
        read_depth,
        &entries,
        inputs.gt21_probabilities,
        inputs.genotype_probabilities,
        inputs.vlen_heads,
        config,
    ))
}

#[allow(clippy::too_many_arguments)]
fn finish(
    chromosome: &str,
    position: u64,
    reference_sequence: &str,
    tensor_position_center: usize,
    info_flag: &str,
    read_depth: u32,
    entries: &AltInfoEntries,
    gt21_probabilities: &[f32],
    genotype_probabilities: &[f32],
    vlen_heads: Option<(&[f32], &[f32])>,
    config: &OutputWithConfig,
) -> Option<String> {
    let decision = if let Some((vlen1, vlen2)) = vlen_heads {
        output_from_fa(
            reference_sequence,
            tensor_position_center,
            gt21_probabilities,
            genotype_probabilities,
            vlen1,
            vlen2,
            entries,
            config.maximum_variant_length_that_need_infer,
        )
    } else {
        output_from(
            reference_sequence,
            tensor_position_center,
            gt21_probabilities,
            genotype_probabilities,
            entries,
            config.maximum_variant_length_that_need_infer,
        )
    };

    if (!config.is_show_reference && decision.is_reference)
        || (!decision.is_reference && decision.reference_base == decision.alternate_base)
    {
        return None;
    }
    if decision.reference_base.is_empty() || decision.alternate_base.is_empty() {
        return None;
    }

    let is_multi = decision.alternate_base.contains(',');

    if config.is_haploid_precise_mode_enabled
        && (decision.is_hetero_snp
            || decision.is_hetero_acgt_ins
            || decision.is_hetero_insins
            || decision.is_hetero_acgt_del
            || decision.is_hetero_deldel
            || decision.is_insertion_and_deletion)
    {
        return None;
    }
    if config.is_haploid_sensitive_mode_enabled && is_multi {
        return None;
    }

    let mut genotype_string = if decision.is_reference {
        "0/0".to_string()
    } else if decision.is_homo_snp || decision.is_homo_insertion || decision.is_homo_deletion {
        "1/1".to_string()
    } else {

        "0/1".to_string()
    };
    if is_multi {
        genotype_string = "1/2".to_string();
    }

    let (alt_type_split, ref_count) = split_alt_info(entries);
    let (supported_reads_count, alt_list_count) = compute_ad(
        &decision,
        &decision.alternate_base,
        &decision.reference_base,
        &alt_type_split,
    );

    let mut alternate_for_emit = decision.alternate_base.clone();
    let mut alt_list_count_out = alt_list_count.clone();
    let mut supported_for_emit = supported_reads_count;

    if decision.is_reference {
        supported_for_emit = ref_count;
        alternate_for_emit = ".".to_string();
        alt_list_count_out.clear();
    }

    let mut allele_frequency = if read_depth > 0 {
        f64::from(supported_for_emit) / f64::from(read_depth)
    } else {
        0.0
    };
    if allele_frequency > 1.0 {
        allele_frequency = 1.0;
    }

    let quality_score = quality_score_from(decision.maximum_probability);
    let filtration_value = filtration_value_from(
        config.quality_score_for_pass,
        quality_score,
        decision.is_reference,
    );

    let (mut reference_base, mut alternate_base) = if !config.keep_iupac_bases {
        (
            convert_iupac_to_n(&decision.reference_base),
            convert_iupac_to_n(&alternate_for_emit),
        )
    } else {
        (decision.reference_base.clone(), alternate_for_emit)
    };
    if decision.is_reference {
        alternate_base = ".".to_string();
    }
    if reference_base.is_empty() {

        reference_base = decision.reference_base.clone();
    }

    let allele_depth = if alt_list_count_out.is_empty() {
        ref_count.to_string()
    } else {
        let mut s = ref_count.to_string();
        for c in &alt_list_count_out {
            s.push(',');
            s.push_str(&c.to_string());
        }
        s
    };
    let allele_frequency_s = render_af(&alt_list_count_out, read_depth, allele_frequency);

    let row = format!(
        "{chromosome}\t{position}\t.\t{reference_base}\t{alternate_base}\t{quality_score:.2}\t{filtration_value}\t{info_flag}\tGT:GQ:DP:AD:AF\t{genotype_string}:{quality_int}:{read_depth}:{allele_depth}:{allele_frequency_s}\n",
        quality_int = quality_score as i64,
    );
    Some(row)
}

#[derive(Debug, Default)]
struct AltTypeSplit {
    snp_map: std::collections::HashMap<String, u32>,
    ins_entries: Vec<(String, u32)>,
    del_entries: Vec<(String, u32)>,
}

impl AltTypeSplit {
    fn ins_count(&self, key: &str) -> u32 {
        self.ins_entries
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, c)| *c)
            .unwrap_or(0)
    }
    fn del_count(&self, key: &str) -> u32 {
        self.del_entries
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, c)| *c)
            .unwrap_or(0)
    }

    fn first_del_with_len(&self, len: usize) -> u32 {
        self.del_entries
            .iter()
            .find(|(k, _)| k.len() == len)
            .map(|(_, c)| *c)
            .unwrap_or(0)
    }
}

fn split_alt_info(entries: &AltInfoEntries) -> (AltTypeSplit, u32) {
    let mut out = AltTypeSplit::default();
    let mut ref_count: u32 = 0;
    for (k, v) in entries {
        let bytes = k.as_bytes();
        if bytes.is_empty() {
            continue;
        }
        match bytes[0] {
            b'X' if bytes.len() >= 2 => {
                out.snp_map.insert((bytes[1] as char).to_string(), *v);
            }
            b'I' => {
                out.ins_entries.push((k[1..].to_string(), *v));
            }
            b'D' => {
                out.del_entries.push((k[1..].to_string(), *v));
            }
            b'R' => {
                ref_count = *v;
            }
            _ => {}
        }
    }
    (out, ref_count)
}

fn compute_ad(
    d: &OutputDecision,
    alternate_base: &str,
    reference_base: &str,
    alt_type_split: &AltTypeSplit,
) -> (u32, Vec<u32>) {
    let mut supported: u32 = 0;
    let mut counts: Vec<u32> = Vec::new();

    if d.is_homo_snp || d.is_hetero_snp {
        for base in alternate_base.chars() {
            if base == ',' {
                continue;
            }
            let c = *alt_type_split.snp_map.get(&base.to_string()).unwrap_or(&0);
            counts.push(c);
            supported += c;
        }
    } else if d.is_homo_insertion || d.is_hetero_insins {
        for ins in alternate_base.split(',') {
            let c = alt_type_split.ins_count(ins);
            counts.push(c);
            supported += c;
        }
    } else if d.is_hetero_acgt_ins {
        let is_snp_ins_multi = alternate_base.contains(',');
        let (snp_base, ins_bases) = if is_snp_ins_multi {
            let parts: Vec<&str> = alternate_base.splitn(2, ',').collect();

            let snp_base = parts[0].chars().next().map(|c| c.to_string());
            let ins_bases = parts[1].to_string();
            (snp_base, ins_bases)
        } else {
            (None, alternate_base.to_string())
        };
        let snp_count = snp_base
            .as_ref()
            .and_then(|b| alt_type_split.snp_map.get(b))
            .copied()
            .unwrap_or(0);
        let ins_count = alt_type_split.ins_count(&ins_bases);
        supported = ins_count + snp_count;
        if snp_base.is_some() {
            counts.push(snp_count);
        }
        counts.push(ins_count);
    } else if d.is_homo_deletion {
        let del_bases = if reference_base.len() > 1 {
            &reference_base[1..]
        } else {
            ""
        };
        let c = alt_type_split.del_count(del_bases);
        counts.push(c);
        supported = c;
    } else if d.is_hetero_deldel {
        if alt_type_split.del_entries.len() > 1 {
            for bases in alternate_base.split(',') {

                let alt_len = reference_base.len().saturating_sub(bases.len());
                let count = alt_type_split.first_del_with_len(alt_len);
                counts.push(count);
                supported += count;
            }
        }
    } else if d.is_hetero_acgt_del {
        let alt_list: Vec<&str> = alternate_base.split(',').collect();
        let is_snp_del_multi = !alt_list.is_empty() && alternate_base.contains(',');
        let snp_base = if is_snp_del_multi && alt_list.len() > 1 {
            alt_list[1].chars().next().map(|c| c.to_string())
        } else {
            None
        };
        let del_bases = if reference_base.len() > 1 {
            &reference_base[1..]
        } else {
            ""
        };
        let snp_count = snp_base
            .as_ref()
            .and_then(|b| alt_type_split.snp_map.get(b))
            .copied()
            .unwrap_or(0);
        let del_count = alt_type_split.del_count(del_bases);
        supported = del_count + snp_count;
        if snp_base.is_some() {
            counts.push(snp_count);
        }
        counts.push(del_count);
    } else if d.is_insertion_and_deletion {
        for bases in alternate_base.split(',') {
            let alt_len_signed = reference_base.len() as i64 - bases.len() as i64;
            let count = if alt_len_signed < 0 {

                let ins_bases = if reference_base.len() > 1 {
                    let cut = bases.len().saturating_sub(reference_base.len() - 1);
                    &bases[..cut]
                } else {
                    bases
                };
                alt_type_split.ins_count(ins_bases)
            } else {
                alt_type_split.first_del_with_len(alt_len_signed as usize)
            };
            counts.push(count);
            supported += count;
        }
    }
    (supported, counts)
}

fn render_af(alt_list_count: &[u32], read_depth: u32, allele_frequency: f64) -> String {
    if alt_list_count.len() <= 1 {
        format!("{:.4}", allele_frequency)
    } else {
        let mut s = String::new();
        for (i, c) in alt_list_count.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let af = if read_depth > 0 {
                (f64::from(*c) / f64::from(read_depth)).min(1.0)
            } else {
                0.0
            };
            s.push_str(&format!("{:.4}", af));
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg() -> OutputWithConfig {
        OutputWithConfig::default()
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
    fn refcall_filtered_when_show_ref_false() {
        let gt21 = one_hot_gt21(0, 0.95);
        let geno = one_hot_geno(0, 0.9);
        let res = output_with(
            OutputWithInputs {
                chr_pos_seq: "chr20:100:A",
                alt_info: "33-RA 33",
                gt21_probabilities: &gt21,
                genotype_probabilities: &geno,
            vlen_heads: None,
        },
            &cfg(),
        )
        .unwrap();
        assert!(res.is_none());
    }

    #[test]
    fn refcall_emitted_when_show_ref_true() {
        let gt21 = one_hot_gt21(0, 0.95);
        let geno = one_hot_geno(0, 0.9);
        let mut c = cfg();
        c.is_show_reference = true;
        let res = output_with(
            OutputWithInputs {
                chr_pos_seq: "chr20:100:A",
                alt_info: "33-RA 33",
                gt21_probabilities: &gt21,
                genotype_probabilities: &geno,
            vlen_heads: None,
        },
            &c,
        )
        .unwrap();
        let row = res.expect("expected ref row");

        assert!(row.starts_with("chr20\t100\t."));
        assert!(row.contains("\tRefCall\tP\t"));
        assert!(row.contains("\t0/0:"));
        assert!(row.ends_with('\n'));
    }

    #[test]
    fn homo_snp_emits_correct_row() {
        use super::super::labels::Gt21;
        let gt21 = one_hot_gt21(Gt21::CC as usize, 0.9);
        let geno = one_hot_geno(1, 0.95);
        let res = output_with(
            OutputWithInputs {
                chr_pos_seq: "chr20:100:A",
                alt_info: "30-XC 30",
                gt21_probabilities: &gt21,
                genotype_probabilities: &geno,
            vlen_heads: None,
        },
            &cfg(),
        )
        .unwrap()
        .unwrap();

        assert!(res.starts_with("chr20\t100\t.\tA\tC\t"), "got: {res}");
        assert!(res.contains("\tPASS\tP\tGT:GQ:DP:AD:AF\t1/1:"), "got: {res}");
        assert!(res.contains(":30:0,30:1.0000\n"));
    }

    #[test]
    fn hetero_snp_emits_correct_row() {
        use super::super::labels::Gt21;
        let gt21 = one_hot_gt21(Gt21::AC as usize, 0.85);
        let geno = one_hot_geno(2, 0.9);
        let res = output_with(
            OutputWithInputs {
                chr_pos_seq: "chr20:100:A",
                alt_info: "33-XC 16 RA 17",
                gt21_probabilities: &gt21,
                genotype_probabilities: &geno,
            vlen_heads: None,
        },
            &cfg(),
        )
        .unwrap()
        .unwrap();
        assert!(res.contains("\tA\tC\t"));
        assert!(res.contains("\t0/1:"));

        assert!(res.contains(":17,16:"));
    }
}
