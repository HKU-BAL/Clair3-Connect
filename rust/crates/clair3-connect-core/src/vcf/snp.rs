use super::labels::{HETERO_SNP_LABELS, HOMO_SNP_LABELS};

pub fn homo_snp_bases_from(homo_snp_probabilities: &[f32]) -> (char, char) {
    assert_eq!(
        homo_snp_probabilities.len(),
        HOMO_SNP_LABELS.len(),
        "homo_snp_probabilities length"
    );
    let idx = argmax(homo_snp_probabilities);
    let label = HOMO_SNP_LABELS[idx];
    let bytes = label.as_bytes();
    (bytes[0] as char, bytes[1] as char)
}

pub fn hetero_snp_bases_from(hetero_snp_probabilities: &[f32]) -> (char, char) {
    assert_eq!(
        hetero_snp_probabilities.len(),
        HETERO_SNP_LABELS.len(),
        "hetero_snp_probabilities length"
    );
    let idx = argmax(hetero_snp_probabilities);
    let label = HETERO_SNP_LABELS[idx];
    let bytes = label.as_bytes();
    (bytes[0] as char, bytes[1] as char)
}

fn argmax(v: &[f32]) -> usize {
    for (i, &x) in v.iter().enumerate() {
        if x.is_nan() {
            return i;
        }
    }
    let mut best_i = 0;
    let mut best_v = f32::NEG_INFINITY;
    for (i, &x) in v.iter().enumerate() {
        if x > best_v {
            best_v = x;
            best_i = i;
        }
    }
    best_i
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn homo_snp_picks_aa() {
        let p = [0.9, 0.05, 0.03, 0.02];
        assert_eq!(homo_snp_bases_from(&p), ('A', 'A'));
    }

    #[test]
    fn homo_snp_picks_tt() {
        let p = [0.01, 0.02, 0.07, 0.9];
        assert_eq!(homo_snp_bases_from(&p), ('T', 'T'));
    }

    #[test]
    fn hetero_snp_picks_ag() {

        let p = [0.1, 0.7, 0.05, 0.05, 0.05, 0.05];
        assert_eq!(hetero_snp_bases_from(&p), ('A', 'G'));
    }

    #[test]
    fn hetero_snp_picks_ct() {
        let p = [0.05, 0.05, 0.05, 0.1, 0.6, 0.15];
        assert_eq!(hetero_snp_bases_from(&p), ('C', 'T'));
    }

    #[test]
    fn argmax_returns_first_max_on_tie() {
        assert_eq!(argmax(&[0.5, 0.5, 0.1]), 0);
        assert_eq!(argmax(&[0.1, 0.5, 0.5]), 1);
    }

    #[test]
    fn argmax_first_nan_wins() {

        assert_eq!(argmax(&[f32::NAN, 0.3, 0.2]), 0);
        assert_eq!(argmax(&[0.5, f32::NAN, 0.9]), 1);

        assert_eq!(argmax(&[0.5, 0.9, 0.9]), 1);
    }

    #[test]
    #[should_panic(expected = "homo_snp_probabilities length")]
    fn homo_snp_panics_on_wrong_length() {
        let _ = homo_snp_bases_from(&[0.1, 0.9]);
    }

    #[test]
    #[should_panic(expected = "hetero_snp_probabilities length")]
    fn hetero_snp_panics_on_wrong_length() {
        let _ = hetero_snp_bases_from(&[0.1, 0.9]);
    }
}
