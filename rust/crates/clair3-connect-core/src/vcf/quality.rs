pub const PHRED_TRANS: f64 = -4.3429448190325175;

pub fn quality_score_from(probability: f32) -> f64 {
    let p = probability as f64;
    let tmp = (PHRED_TRANS * (((1.0 - p) + 1e-10) / (p + 1e-10)).ln() + 10.0).max(0.0);
    (tmp * 100.0).round_ties_even() / 100.0
}

pub fn filtration_value_from(
    quality_score_for_pass: Option<f64>,
    quality_score: f64,
    is_reference: bool,
) -> &'static str {
    if is_reference {
        return "RefCall";
    }
    match quality_score_for_pass {
        None => "PASS",
        Some(threshold) if quality_score >= threshold => "PASS",
        Some(_) => "LowQual",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phred_trans_value() {

        let expected = -10.0_f64 * std::f64::consts::E.log10();
        assert!(
            (PHRED_TRANS - expected).abs() < 1e-14,
            "PHRED_TRANS={PHRED_TRANS} vs expected={expected}"
        );
    }

    #[test]
    fn phred_trans_bit_equal_to_python() {

        let python_value = -10.0_f64 / std::f64::consts::LN_10;
        assert_eq!(
            PHRED_TRANS.to_bits(),
            python_value.to_bits(),
            "PHRED_TRANS literal {PHRED_TRANS} differs from -10/ln(10) = {python_value}"
        );
    }

    #[test]
    fn quality_score_high_probability() {

        let q_high = quality_score_from(0.999_f32);
        let q_mid = quality_score_from(0.9_f32);
        let q_low = quality_score_from(0.6_f32);
        assert!(q_high > q_mid, "q_high={q_high} q_mid={q_mid}");
        assert!(q_mid > q_low, "q_mid={q_mid} q_low={q_low}");
        assert!(q_high > 20.0, "q_high={q_high} should be >20");
    }

    #[test]
    fn quality_score_at_p_half() {

        let q = quality_score_from(0.5_f32);
        assert!((q - 10.0).abs() < 1e-2, "got {q}");
    }

    #[test]
    fn quality_score_clamped_at_zero() {

        let q = quality_score_from(0.0001_f32);
        assert!(q >= 0.0);
    }

    #[test]
    fn filtration_reference_takes_precedence() {
        assert_eq!(
            filtration_value_from(Some(2.0), 100.0, true),
            "RefCall",
            "is_reference overrides quality"
        );
    }

    #[test]
    fn filtration_pass_when_threshold_unset() {
        assert_eq!(filtration_value_from(None, 0.5, false), "PASS");
        assert_eq!(filtration_value_from(None, 100.0, false), "PASS");
    }

    #[test]
    fn filtration_pass_above_threshold() {
        assert_eq!(filtration_value_from(Some(2.0), 5.0, false), "PASS");
        assert_eq!(filtration_value_from(Some(2.0), 2.0, false), "PASS");
    }

    #[test]
    fn filtration_lowqual_below_threshold() {
        assert_eq!(filtration_value_from(Some(2.0), 1.99, false), "LowQual");
    }
}
