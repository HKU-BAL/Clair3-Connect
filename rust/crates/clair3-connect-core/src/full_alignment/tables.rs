pub const NUM2COUNTBASE_FA: [i8; 32] = [
    100, 0, 25, -100, 0, 0, 75, 0,
    -50, 0, 0, 0, 0, 100, 0, 0,
    0, 0, 0, 50, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
];

pub const ACGT2NUM: [i8; 32] = [
    0, 0, 1, 0, 0, 0, 2, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 3, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
];

pub const SEQ_NT16_STR: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

#[inline]
pub fn fa_feature_for_char(c: u8) -> i8 {
    NUM2COUNTBASE_FA[((c.wrapping_sub(b'A')) as usize) & 0x1F]
}

#[inline]
pub fn fa_feature_for_seqi(seqi: u8) -> i8 {
    fa_feature_for_char(SEQ_NT16_STR[(seqi as usize) & 0x0F])
}

#[inline]
pub fn seqi_to_char(seqi: u8) -> u8 {
    SEQ_NT16_STR[(seqi as usize) & 0x0F]
}

#[inline]
pub fn acgt_num_for_char(c: u8) -> i8 {
    ACGT2NUM[((c.wrapping_sub(b'A')) as usize) & 0x1F]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fa_feature_acgtn_letters() {
        assert_eq!(fa_feature_for_char(b'A'), 100);
        assert_eq!(fa_feature_for_char(b'C'), 25);
        assert_eq!(fa_feature_for_char(b'G'), 75);
        assert_eq!(fa_feature_for_char(b'T'), 50);
        assert_eq!(fa_feature_for_char(b'N'), 100);

        assert_eq!(fa_feature_for_char(b'I'), -50);
        assert_eq!(fa_feature_for_char(b'D'), -100);
    }

    #[test]
    fn fa_feature_for_seqi_acgt() {

        assert_eq!(fa_feature_for_seqi(1), 100);
        assert_eq!(fa_feature_for_seqi(2), 25);
        assert_eq!(fa_feature_for_seqi(4), 75);
        assert_eq!(fa_feature_for_seqi(8), 50);
        assert_eq!(fa_feature_for_seqi(15), 100);
    }

    #[test]
    fn seqi_decode_acgt() {
        assert_eq!(seqi_to_char(1), b'A');
        assert_eq!(seqi_to_char(2), b'C');
        assert_eq!(seqi_to_char(4), b'G');
        assert_eq!(seqi_to_char(8), b'T');
        assert_eq!(seqi_to_char(15), b'N');
    }

    #[test]
    fn acgt_num_letters() {
        assert_eq!(acgt_num_for_char(b'A'), 0);
        assert_eq!(acgt_num_for_char(b'C'), 1);
        assert_eq!(acgt_num_for_char(b'G'), 2);
        assert_eq!(acgt_num_for_char(b'T'), 3);
    }
}
