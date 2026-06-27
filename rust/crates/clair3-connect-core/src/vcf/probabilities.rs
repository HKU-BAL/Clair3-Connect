use super::labels::{variant_length, Genotype, Gt21, HETERO_SNP_GT21, HOMO_SNP_GT21};
use super::length_tuples::{
    hetero_del_tuples, hetero_deldel_tuples, hetero_ins_tuples, hetero_insdel_tuples,
    hetero_insins_tuples, homo_del_tuples, homo_ins_tuples,
};

#[derive(Debug, Clone)]
pub enum Outcome {

    HomoRefShortcut(f32),

    Full(Box<FullOutcome>),
}

#[derive(Debug, Clone)]
pub struct FullOutcome {

    pub homo_ref: f32,

    pub homo_snp: [f32; 4],

    pub hetero_snp: [f32; 6],

    pub homo_ins_probabilities: Vec<f32>,

    pub homo_ins_lengths: Vec<i32>,

    pub hetero_insins_probabilities: Vec<f32>,

    pub hetero_insins_length_tuples: Vec<(i32, i32)>,

    pub hetero_acgt_ins_probabilities: Vec<f32>,

    pub hetero_acgt_ins_bases: Vec<char>,

    pub hetero_acgt_ins_lengths: Vec<i32>,

    pub homo_del_probabilities: Vec<f32>,

    pub homo_del_lengths: Vec<i32>,

    pub hetero_deldel_probabilities: Vec<f32>,

    pub hetero_deldel_length_tuples: Vec<(i32, i32)>,

    pub hetero_acgt_del_probabilities: Vec<f32>,

    pub hetero_acgt_del_bases: Vec<char>,

    pub hetero_acgt_del_lengths: Vec<i32>,

    pub hetero_insdel_probabilities: Vec<f32>,

    pub hetero_insdel_length_tuples: Vec<(i32, i32)>,
}

pub fn possible_outcomes_pileup(
    gt21_probabilities: &[f32],
    genotype_probabilities: &[f32],
    reference_base: u8,
) -> Outcome {
    debug_assert_eq!(gt21_probabilities.len(), 21);
    debug_assert_eq!(genotype_probabilities.len(), 3);

    let homo_ref_p = genotype_probabilities[Genotype::HomoReference as usize];
    let homo_var_p = genotype_probabilities[Genotype::HomoVariant as usize];
    let het_var_p = genotype_probabilities[Genotype::HeteroVariant as usize];

    let ref_label_owned = format!("{}{}", reference_base as char, reference_base as char);
    let reference_gt21 = Gt21::from_label(&ref_label_owned)
        .expect("reference_base must be ACGT after iupac_to_acgt");
    let ref_gt21_p = gt21_probabilities[reference_gt21 as usize];

    let homo_ref_probability = homo_ref_p * ref_gt21_p;

    if homo_ref_p >= 0.5 && ref_gt21_p >= 0.5 {
        return Outcome::HomoRefShortcut(homo_ref_probability);
    }

    let mut homo_snp = [0.0_f32; 4];
    for (i, gt21) in HOMO_SNP_GT21.iter().enumerate() {
        homo_snp[i] = homo_var_p * gt21_probabilities[*gt21 as usize];
    }
    let mut hetero_snp = [0.0_f32; 6];
    for (i, gt21) in HETERO_SNP_GT21.iter().enumerate() {
        hetero_snp[i] = het_var_p * gt21_probabilities[*gt21 as usize];
    }

    let homo_ins_probabilities = vec![homo_var_p * gt21_probabilities[Gt21::InsIns as usize]];
    let hetero_insins_probabilities =
        vec![het_var_p * gt21_probabilities[Gt21::InsIns as usize]];
    let hetero_acgt_ins_probabilities: Vec<f32> = [
        Gt21::AIns,
        Gt21::CIns,
        Gt21::GIns,
        Gt21::TIns,
    ]
    .iter()
    .map(|g| het_var_p * gt21_probabilities[*g as usize])
    .collect();

    let homo_del_probabilities = vec![homo_var_p * gt21_probabilities[Gt21::DelDel as usize]];
    let hetero_deldel_probabilities =
        vec![het_var_p * gt21_probabilities[Gt21::DelDel as usize]];
    let hetero_acgt_del_probabilities: Vec<f32> = [
        Gt21::ADel,
        Gt21::CDel,
        Gt21::GDel,
        Gt21::TDel,
    ]
    .iter()
    .map(|g| het_var_p * gt21_probabilities[*g as usize])
    .collect();

    let hetero_insdel_probabilities = vec![het_var_p * gt21_probabilities[Gt21::InsDel as usize]];

    Outcome::Full(Box::new(FullOutcome {
        homo_ref: homo_ref_probability,
        homo_snp,
        hetero_snp,
        homo_ins_probabilities,
        homo_ins_lengths: Vec::new(),
        hetero_insins_probabilities,
        hetero_insins_length_tuples: Vec::new(),
        hetero_acgt_ins_probabilities,
        hetero_acgt_ins_bases: Vec::new(),
        hetero_acgt_ins_lengths: Vec::new(),
        homo_del_probabilities,
        homo_del_lengths: Vec::new(),
        hetero_deldel_probabilities,
        hetero_deldel_length_tuples: Vec::new(),
        hetero_acgt_del_probabilities,
        hetero_acgt_del_bases: Vec::new(),
        hetero_acgt_del_lengths: Vec::new(),
        hetero_insdel_probabilities,
        hetero_insdel_length_tuples: Vec::new(),
    }))
}

pub fn possible_outcomes_fa(
    gt21_probabilities: &[f32],
    genotype_probabilities: &[f32],
    vlen1: &[f32],
    vlen2: &[f32],
    reference_base: u8,
) -> Outcome {
    debug_assert_eq!(gt21_probabilities.len(), 21);
    debug_assert_eq!(genotype_probabilities.len(), 3);
    debug_assert_eq!(vlen1.len(), variant_length::OUTPUT_LABEL_COUNT);
    debug_assert_eq!(vlen2.len(), variant_length::OUTPUT_LABEL_COUNT);

    let homo_ref_p = genotype_probabilities[Genotype::HomoReference as usize];
    let homo_var_p = genotype_probabilities[Genotype::HomoVariant as usize];
    let het_var_p = genotype_probabilities[Genotype::HeteroVariant as usize];

    let ref_label_owned = format!("{}{}", reference_base as char, reference_base as char);
    let reference_gt21 = Gt21::from_label(&ref_label_owned)
        .expect("reference_base must be ACGT after iupac_to_acgt");
    let ref_gt21_p = gt21_probabilities[reference_gt21 as usize];

    let offset = variant_length::INDEX_OFFSET;
    let vl1_0 = vlen1[offset];
    let vl2_0 = vlen2[offset];
    let vlen_0_p = vl1_0 * vl2_0;
    let homo_ref_probability = vlen_0_p * homo_ref_p * ref_gt21_p;

    if vl1_0 >= 0.5 && vl2_0 >= 0.5 && homo_ref_p >= 0.5 && ref_gt21_p >= 0.5 {
        return Outcome::HomoRefShortcut(homo_ref_probability);
    }

    let mut homo_snp = [0.0_f32; 4];
    for (i, gt21) in HOMO_SNP_GT21.iter().enumerate() {
        homo_snp[i] = vlen_0_p * homo_var_p * gt21_probabilities[*gt21 as usize];
    }
    let mut hetero_snp = [0.0_f32; 6];
    for (i, gt21) in HETERO_SNP_GT21.iter().enumerate() {
        hetero_snp[i] = vlen_0_p * het_var_p * gt21_probabilities[*gt21 as usize];
    }

    let homo_ins_extra = homo_var_p * gt21_probabilities[Gt21::InsIns as usize];
    let homo_ins = homo_ins_tuples(vlen1, vlen2, homo_ins_extra);
    let homo_ins_lengths: Vec<i32> = homo_ins.iter().map(|(l, _)| *l).collect();
    let homo_ins_probabilities: Vec<f32> = homo_ins.iter().map(|(_, p)| *p).collect();

    let hetero_insins_extra = het_var_p * gt21_probabilities[Gt21::InsIns as usize];
    let hetero_insins = hetero_insins_tuples(vlen1, vlen2, hetero_insins_extra);
    let hetero_insins_length_tuples: Vec<(i32, i32)> = hetero_insins.iter().map(|(t, _)| *t).collect();
    let hetero_insins_probabilities: Vec<f32> = hetero_insins.iter().map(|(_, p)| *p).collect();

    let hetero_ins_base_lens = hetero_ins_tuples(vlen1, vlen2);
    let mut hetero_acgt_ins_bases: Vec<char> = Vec::new();
    let mut hetero_acgt_ins_lengths: Vec<i32> = Vec::new();
    let mut hetero_acgt_ins_probabilities: Vec<f32> = Vec::new();
    let acgt_ins_pairs: [(Gt21, char); 4] = [
        (Gt21::AIns, 'A'),
        (Gt21::CIns, 'C'),
        (Gt21::GIns, 'G'),
        (Gt21::TIns, 'T'),
    ];
    for (length, len_p) in &hetero_ins_base_lens {
        for (gt21, base_char) in &acgt_ins_pairs {
            hetero_acgt_ins_bases.push(*base_char);
            hetero_acgt_ins_lengths.push(*length);
            hetero_acgt_ins_probabilities
                .push(len_p * gt21_probabilities[*gt21 as usize] * het_var_p);
        }
    }

    let homo_del_extra = homo_var_p * gt21_probabilities[Gt21::DelDel as usize];
    let homo_del = homo_del_tuples(vlen1, vlen2, homo_del_extra);
    let homo_del_lengths: Vec<i32> = homo_del.iter().map(|(l, _)| *l).collect();
    let homo_del_probabilities: Vec<f32> = homo_del.iter().map(|(_, p)| *p).collect();

    let hetero_deldel_extra = het_var_p * gt21_probabilities[Gt21::DelDel as usize];
    let hetero_deldel = hetero_deldel_tuples(vlen1, vlen2, hetero_deldel_extra);
    let hetero_deldel_length_tuples: Vec<(i32, i32)> = hetero_deldel.iter().map(|(t, _)| *t).collect();
    let hetero_deldel_probabilities: Vec<f32> = hetero_deldel.iter().map(|(_, p)| *p).collect();

    let hetero_del_base_lens = hetero_del_tuples(vlen1, vlen2);
    let mut hetero_acgt_del_bases: Vec<char> = Vec::new();
    let mut hetero_acgt_del_lengths: Vec<i32> = Vec::new();
    let mut hetero_acgt_del_probabilities: Vec<f32> = Vec::new();
    let acgt_del_pairs: [(Gt21, char); 4] = [
        (Gt21::ADel, 'A'),
        (Gt21::CDel, 'C'),
        (Gt21::GDel, 'G'),
        (Gt21::TDel, 'T'),
    ];
    for (length, len_p) in &hetero_del_base_lens {
        for (gt21, base_char) in &acgt_del_pairs {
            hetero_acgt_del_bases.push(*base_char);
            hetero_acgt_del_lengths.push(*length);
            hetero_acgt_del_probabilities
                .push(len_p * gt21_probabilities[*gt21 as usize] * het_var_p);
        }
    }

    let hetero_insdel_extra = het_var_p * gt21_probabilities[Gt21::InsDel as usize];
    let hetero_insdel = hetero_insdel_tuples(vlen1, vlen2, hetero_insdel_extra);
    let hetero_insdel_length_tuples: Vec<(i32, i32)> = hetero_insdel.iter().map(|(t, _)| *t).collect();
    let hetero_insdel_probabilities: Vec<f32> = hetero_insdel.iter().map(|(_, p)| *p).collect();

    Outcome::Full(Box::new(FullOutcome {
        homo_ref: homo_ref_probability,
        homo_snp,
        hetero_snp,
        homo_ins_probabilities,
        homo_ins_lengths,
        hetero_insins_probabilities,
        hetero_insins_length_tuples,
        hetero_acgt_ins_probabilities,
        hetero_acgt_ins_bases,
        hetero_acgt_ins_lengths,
        homo_del_probabilities,
        homo_del_lengths,
        hetero_deldel_probabilities,
        hetero_deldel_length_tuples,
        hetero_acgt_del_probabilities,
        hetero_acgt_del_bases,
        hetero_acgt_del_lengths,
        hetero_insdel_probabilities,
        hetero_insdel_length_tuples,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gt21_peak(peak: Gt21, peak_value: f32) -> Vec<f32> {
        let mut v = vec![(1.0 - peak_value) / 20.0; 21];
        v[peak as usize] = peak_value;
        v
    }

    fn geno_peak(peak: Genotype, peak_value: f32) -> Vec<f32> {
        let mut v = vec![(1.0 - peak_value) / 2.0; 3];
        v[peak as usize] = peak_value;
        v
    }

    #[test]
    fn early_exit_when_clearly_homo_ref() {

        let gt21 = gt21_peak(Gt21::AA, 0.95);
        let geno = geno_peak(Genotype::HomoReference, 0.9);
        let out = possible_outcomes_pileup(&gt21, &geno, b'A');
        match out {
            Outcome::HomoRefShortcut(p) => {

                assert!((p - 0.855).abs() < 1e-5, "got {p}");
            }
            _ => panic!("expected HomoRefShortcut"),
        }
    }

    #[test]
    fn full_path_when_homo_ref_low() {

        let gt21 = gt21_peak(Gt21::AC, 0.7);
        let geno = geno_peak(Genotype::HeteroVariant, 0.8);
        let out = possible_outcomes_pileup(&gt21, &geno, b'A');
        match out {
            Outcome::Full(f) => {

                assert!((f.hetero_snp[0] - 0.56).abs() < 1e-5, "got {}", f.hetero_snp[0]);

                assert_eq!(f.homo_ins_probabilities.len(), 1);
                assert_eq!(f.homo_ins_lengths.len(), 0);
                assert_eq!(f.hetero_acgt_ins_probabilities.len(), 4);
                assert_eq!(f.hetero_acgt_ins_bases.len(), 0);
                assert_eq!(f.hetero_insdel_probabilities.len(), 1);
            }
            _ => panic!("expected Full"),
        }
    }

    #[test]
    fn full_path_when_ref_gt21_low() {

        let mut gt21 = vec![0.05_f32; 21];
        gt21[Gt21::AC as usize] = 0.5;
        gt21[Gt21::AA as usize] = 0.45;
        let geno = geno_peak(Genotype::HomoReference, 0.9);
        let out = possible_outcomes_pileup(&gt21, &geno, b'A');
        assert!(matches!(out, Outcome::Full(_)));
    }

    #[test]
    fn homo_snp_ordering_aa_cc_gg_tt() {

        let mut gt21 = vec![0.01_f32; 21];
        gt21[Gt21::AA as usize] = 0.8;
        gt21[Gt21::CC as usize] = 0.04;
        gt21[Gt21::GG as usize] = 0.02;
        gt21[Gt21::TT as usize] = 0.13;
        let geno = geno_peak(Genotype::HomoVariant, 0.9);
        let out = possible_outcomes_pileup(&gt21, &geno, b'C');
        match out {
            Outcome::Full(f) => {

                assert!((f.homo_snp[0] - 0.9 * 0.8).abs() < 1e-5);
                assert!((f.homo_snp[1] - 0.9 * 0.04).abs() < 1e-5);
                assert!((f.homo_snp[2] - 0.9 * 0.02).abs() < 1e-5);
                assert!((f.homo_snp[3] - 0.9 * 0.13).abs() < 1e-5);
            }
            _ => panic!("expected Full"),
        }
    }

    #[test]
    fn hetero_acgt_ins_ordering_a_c_g_t() {
        let mut gt21 = vec![0.01_f32; 21];
        gt21[Gt21::AIns as usize] = 0.1;
        gt21[Gt21::CIns as usize] = 0.2;
        gt21[Gt21::GIns as usize] = 0.3;
        gt21[Gt21::TIns as usize] = 0.05;
        let geno = geno_peak(Genotype::HeteroVariant, 0.7);
        let out = possible_outcomes_pileup(&gt21, &geno, b'A');
        match out {
            Outcome::Full(f) => {
                let v = &f.hetero_acgt_ins_probabilities;
                assert_eq!(v.len(), 4);
                assert!((v[0] - 0.7 * 0.1).abs() < 1e-5);
                assert!((v[1] - 0.7 * 0.2).abs() < 1e-5);
                assert!((v[2] - 0.7 * 0.3).abs() < 1e-5);
                assert!((v[3] - 0.7 * 0.05).abs() < 1e-5);
            }
            _ => panic!("expected Full"),
        }
    }

    #[test]
    fn nan_in_gt21_propagates_no_panic() {

        let mut gt21 = vec![0.01_f32; 21];
        gt21[Gt21::AC as usize] = f32::NAN;
        gt21[Gt21::AA as usize] = 0.4;
        let geno = geno_peak(Genotype::HeteroVariant, 0.7);
        let out = possible_outcomes_pileup(&gt21, &geno, b'A');
        match out {
            Outcome::Full(f) => {
                assert!(f.hetero_snp[0].is_nan(), "got {}", f.hetero_snp[0]);
                assert!(f.homo_ref.is_finite(), "homo_ref={}", f.homo_ref);
            }
            _ => panic!("expected Full"),
        }
    }

    #[test]
    fn shortcut_blocked_when_homo_ref_p_is_nan() {

        let gt21 = gt21_peak(Gt21::AA, 0.95);
        let mut geno = vec![0.05_f32; 3];
        geno[Genotype::HomoReference as usize] = f32::NAN;
        let out = possible_outcomes_pileup(&gt21, &geno, b'A');

        match out {
            Outcome::Full(f) => assert!(f.homo_ref.is_nan()),
            Outcome::HomoRefShortcut(_) => panic!("shortcut on NaN must not trigger"),
        }
    }
}
