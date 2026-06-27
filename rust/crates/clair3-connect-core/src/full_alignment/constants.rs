pub const HAP_UNPHASED: u8 = 0;

pub const HAP_1: u8 = 1;

pub const HAP_2: u8 = 2;

pub const HAP_TYPE: [i8; 3] = [60, 30, 90];

pub const OVERHANG: usize = 10;

pub const MIN_HAPLOTAG_MQ: u8 = 20;

pub const EXPAND_REFERENCE_REGION: i64 = 2_000_000;

pub const FLANKING_BASE_NUM: usize = 16;

pub const NO_OF_POSITIONS: usize = 33;

pub const CHANNEL_SIZE: usize = 8;

pub const CHANNEL_SIZE_DWELL: usize = 9;

pub const MIN_BQ: u8 = 0;

pub const SAMTOOLS_VIEW_FILTER_FLAG: u16 = 2316;

pub const MAX_READ_COUNT: usize = 1000;

pub const MAX_INDEL_LENGTH: usize = 50;

#[inline]
pub fn normalize_mq(x: u8) -> i8 {
    if x < 60 {
        (100.0_f64 * (x as f64) / 60.0) as i8
    } else {
        100
    }
}

#[inline]
pub fn normalize_bq(x: u8) -> i8 {
    if x < 40 {
        (100.0_f64 * (x as f64) / 40.0) as i8
    } else {
        100
    }
}

#[inline]
pub fn normalize_af(x: f64) -> i8 {
    if x < 1.0 {
        (100.0 * x) as i8
    } else {
        100
    }
}

#[inline]
pub fn normalize_strand(is_reverse: bool) -> i8 {
    if is_reverse {
        50
    } else {
        100
    }
}

#[inline]
pub fn normalize_hap(hap: u8) -> i8 {
    HAP_TYPE[hap as usize]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mq_normalize_examples() {

        assert_eq!(normalize_mq(60), 100);
        assert_eq!(normalize_mq(30), 50);
        assert_eq!(normalize_mq(45), 75);
        assert_eq!(normalize_mq(0), 0);

        assert_eq!(normalize_mq(59), 98);
    }

    #[test]
    fn bq_normalize_examples() {
        assert_eq!(normalize_bq(40), 100);
        assert_eq!(normalize_bq(20), 50);
        assert_eq!(normalize_bq(39), 97);
    }

    #[test]
    fn af_normalize_examples() {
        assert_eq!(normalize_af(1.0), 100);
        assert_eq!(normalize_af(0.5), 50);
        assert_eq!(normalize_af(0.99), 99);
        assert_eq!(normalize_af(0.0), 0);
    }

    #[test]
    fn strand_normalize() {
        assert_eq!(normalize_strand(false), 100);
        assert_eq!(normalize_strand(true), 50);
    }

    #[test]
    fn hap_normalize() {
        assert_eq!(normalize_hap(HAP_UNPHASED), 60);
        assert_eq!(normalize_hap(HAP_1), 30);
        assert_eq!(normalize_hap(HAP_2), 90);
    }
}
