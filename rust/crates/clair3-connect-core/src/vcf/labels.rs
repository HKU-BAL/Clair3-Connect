pub const GT21_COUNT: usize = 21;

pub const GT21_LABELS: [&str; GT21_COUNT] = [
    "AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT", "DelDel", "ADel", "CDel", "GDel",
    "TDel", "InsIns", "AIns", "CIns", "GIns", "TIns", "InsDel",
];

#[allow(non_camel_case_types, missing_docs)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Gt21 {
    AA = 0,
    AC = 1,
    AG = 2,
    AT = 3,
    CC = 4,
    CG = 5,
    CT = 6,
    GG = 7,
    GT = 8,
    TT = 9,
    DelDel = 10,
    ADel = 11,
    CDel = 12,
    GDel = 13,
    TDel = 14,
    InsIns = 15,
    AIns = 16,
    CIns = 17,
    GIns = 18,
    TIns = 19,
    InsDel = 20,
}

impl Gt21 {

    pub fn from_usize(i: usize) -> Option<Self> {
        Some(match i {
            0 => Self::AA,
            1 => Self::AC,
            2 => Self::AG,
            3 => Self::AT,
            4 => Self::CC,
            5 => Self::CG,
            6 => Self::CT,
            7 => Self::GG,
            8 => Self::GT,
            9 => Self::TT,
            10 => Self::DelDel,
            11 => Self::ADel,
            12 => Self::CDel,
            13 => Self::GDel,
            14 => Self::TDel,
            15 => Self::InsIns,
            16 => Self::AIns,
            17 => Self::CIns,
            18 => Self::GIns,
            19 => Self::TIns,
            20 => Self::InsDel,
            _ => return None,
        })
    }

    pub fn label(self) -> &'static str {
        GT21_LABELS[self as usize]
    }

    pub fn from_label(s: &str) -> Option<Self> {
        Self::from_usize(GT21_LABELS.iter().position(|&l| l == s)?)
    }
}

pub const HOMO_SNP_GT21: [Gt21; 4] = [Gt21::AA, Gt21::CC, Gt21::GG, Gt21::TT];

pub const HOMO_SNP_LABELS: [&str; 4] = ["AA", "CC", "GG", "TT"];

pub const HETERO_SNP_GT21: [Gt21; 6] = [Gt21::AC, Gt21::AG, Gt21::AT, Gt21::CG, Gt21::CT, Gt21::GT];

pub const HETERO_SNP_LABELS: [&str; 6] = ["AC", "AG", "AT", "CG", "CT", "GT"];

pub fn partial_label_from(reference: &str, alternate: &str) -> String {
    if reference.len() > alternate.len() {
        "Del".to_string()
    } else if reference.len() < alternate.len() {
        "Ins".to_string()
    } else {

        alternate
            .chars()
            .next()
            .map(|c| c.to_string())
            .unwrap_or_default()
    }
}

pub fn mix_two_partial_labels(label1: &str, label2: &str) -> String {

    if label1.len() == 1 && label2.len() == 1 {
        return if label1 <= label2 {
            format!("{label1}{label2}")
        } else {
            format!("{label2}{label1}")
        };
    }

    if label1.len() > 1 && label2.len() == 1 {
        return format!("{label2}{label1}");
    }
    if label2.len() > 1 && label1.len() == 1 {
        return format!("{label1}{label2}");
    }

    if !label1.is_empty() && label1 == label2 {
        return format!("{label1}{label2}");
    }

    "InsDel".to_string()
}

pub const GENOTYPE_STRINGS: [&str; 4] = ["0/0", "1/1", "0/1", "1/2"];

#[allow(missing_docs)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Genotype {

    HomoReference = 0,

    HomoVariant = 1,

    HeteroVariant = 2,

    HeteroVariantMulti = 3,
}

impl Genotype {

    pub fn as_str(self) -> &'static str {
        GENOTYPE_STRINGS[self as usize]
    }

    pub fn from_pair(g1: u32, g2: u32) -> Self {
        if g1 == 0 && g2 == 0 {
            Self::HomoReference
        } else if g1 == g2 {
            Self::HomoVariant
        } else if g1 != 0 && g2 != 0 {
            Self::HeteroVariantMulti
        } else {
            Self::HeteroVariant
        }
    }

    pub fn for_task(self) -> Self {
        if matches!(self, Self::HeteroVariantMulti) {
            Self::HeteroVariant
        } else {
            self
        }
    }
}

pub mod variant_length {

    pub const INDEX_OFFSET: usize = 16;

    pub const MIN: i32 = -(INDEX_OFFSET as i32);

    pub const MAX: i32 = INDEX_OFFSET as i32;

    pub const OUTPUT_LABEL_COUNT: usize = 2 * INDEX_OFFSET + 1;
}

pub const LABEL_SHAPE_CUM: [usize; 4] = [21, 24, 57, 90];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gt21_round_trip() {
        for (i, label) in GT21_LABELS.iter().enumerate() {
            let g = Gt21::from_usize(i).unwrap();
            assert_eq!(g as usize, i);
            assert_eq!(g.label(), *label);
            assert_eq!(Gt21::from_label(label), Some(g));
        }
        assert!(Gt21::from_usize(21).is_none());
        assert!(Gt21::from_label("ZZ").is_none());
    }

    #[test]
    fn homo_hetero_snp_indices() {
        assert_eq!(HOMO_SNP_GT21.len(), 4);
        assert_eq!(HETERO_SNP_GT21.len(), 6);
        for (i, g) in HOMO_SNP_GT21.iter().enumerate() {
            assert_eq!(g.label(), HOMO_SNP_LABELS[i]);
        }
        for (i, g) in HETERO_SNP_GT21.iter().enumerate() {
            assert_eq!(g.label(), HETERO_SNP_LABELS[i]);
        }
    }

    #[test]
    fn partial_label_basics() {
        assert_eq!(partial_label_from("A", "C"), "C");
        assert_eq!(partial_label_from("A", "AG"), "Ins");
        assert_eq!(partial_label_from("AC", "A"), "Del");
        assert_eq!(partial_label_from("A", "A"), "A");
    }

    #[test]
    fn mix_labels_snp_pairs() {

        assert_eq!(mix_two_partial_labels("A", "G"), "AG");
        assert_eq!(mix_two_partial_labels("G", "A"), "AG");
        assert_eq!(mix_two_partial_labels("C", "C"), "CC");
    }

    #[test]
    fn mix_labels_indels() {

        assert_eq!(mix_two_partial_labels("A", "Ins"), "AIns");
        assert_eq!(mix_two_partial_labels("Del", "G"), "GDel");

        assert_eq!(mix_two_partial_labels("Ins", "Ins"), "InsIns");
        assert_eq!(mix_two_partial_labels("Del", "Del"), "DelDel");

        assert_eq!(mix_two_partial_labels("Ins", "Del"), "InsDel");
    }

    #[test]
    fn genotype_pair_mapping() {
        assert_eq!(Genotype::from_pair(0, 0), Genotype::HomoReference);
        assert_eq!(Genotype::from_pair(1, 1), Genotype::HomoVariant);
        assert_eq!(Genotype::from_pair(0, 1), Genotype::HeteroVariant);
        assert_eq!(Genotype::from_pair(1, 0), Genotype::HeteroVariant);
        assert_eq!(Genotype::from_pair(1, 2), Genotype::HeteroVariantMulti);
        assert_eq!(
            Genotype::from_pair(1, 2).for_task(),
            Genotype::HeteroVariant
        );
        assert_eq!(Genotype::HomoVariant.for_task(), Genotype::HomoVariant);
        assert_eq!(Genotype::HomoReference.as_str(), "0/0");
        assert_eq!(Genotype::HomoVariant.as_str(), "1/1");
        assert_eq!(Genotype::HeteroVariant.as_str(), "0/1");
        assert_eq!(Genotype::HeteroVariantMulti.as_str(), "1/2");
    }

    #[test]
    fn variant_length_constants() {
        use variant_length::*;
        assert_eq!(INDEX_OFFSET, 16);
        assert_eq!(MIN, -16);
        assert_eq!(MAX, 16);
        assert_eq!(OUTPUT_LABEL_COUNT, 33);
    }

    #[test]
    fn label_shape_cum_widths() {
        assert_eq!(LABEL_SHAPE_CUM, [21, 24, 57, 90]);
    }
}
