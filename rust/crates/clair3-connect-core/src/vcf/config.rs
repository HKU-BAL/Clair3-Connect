#[derive(Debug, Clone)]
pub struct OutputConfig {

    pub is_show_reference: bool,

    pub is_haploid_precise: bool,

    pub is_haploid_sensitive: bool,

    pub quality_score_for_pass: Option<f64>,

    pub add_indel_length: bool,

    pub gvcf: bool,

    pub pileup: bool,

    pub enable_long_indel: bool,

    pub maximum_variant_length_that_need_infer: usize,

    pub keep_iupac_bases: bool,

    pub is_debug: bool,
}

impl Default for OutputConfig {
    fn default() -> Self {
        Self {
            is_show_reference: true,
            is_haploid_precise: false,
            is_haploid_sensitive: false,
            quality_score_for_pass: Some(2.0),
            add_indel_length: false,
            gvcf: false,
            pileup: true,
            enable_long_indel: false,
            maximum_variant_length_that_need_infer: 50,
            keep_iupac_bases: false,
            is_debug: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_mirrors_stock_pileup_cli() {
        let c = OutputConfig::default();
        assert!(c.is_show_reference);
        assert!(c.pileup);
        assert!(!c.add_indel_length);
        assert!(!c.gvcf);
        assert!(!c.keep_iupac_bases);
        assert_eq!(c.quality_score_for_pass, Some(2.0));
        assert_eq!(c.maximum_variant_length_that_need_infer, 50);
    }
}
