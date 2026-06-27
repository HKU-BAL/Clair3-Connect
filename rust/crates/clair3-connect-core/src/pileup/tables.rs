pub const PLP_BASES_CLAIR3: &[u8; 18] = b"ACGTIIDDDacgtiiddd";

pub const NUM_TO_COUNT_BASE_CLAIR3: [i32; 32] = [
    -1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1,
    -1, 9, 10, -1, 11, -1, -1, -1, 12, -1, -1, -1, -1, -1, -1, -1,
];

pub const BASE2INDEX: [i32; 32] = [
    0, 0, 1, 0, 0, 0, 2, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 3, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
];

pub const SEQ_NT16_STR: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

pub const FEAT_LEN: usize = 18;

pub const C3_FWD_INS_ALL: usize = 4;

pub const C3_FWD_INS_BEST: usize = 5;

pub const C3_FWD_DEL_ALL: usize = 6;

pub const C3_FWD_DEL_BEST: usize = 7;

pub const C3_FWD_DEL: usize = 8;

pub const C3_REV_INS_ALL: usize = 13;

pub const C3_REV_INS_BEST: usize = 14;

pub const C3_REV_DEL_ALL: usize = 15;

pub const C3_REV_DEL_BEST: usize = 16;

pub const C3_REV_DEL: usize = 17;

pub const REVERSE_POS_START: usize = 9;

pub const FLANKING_BASE_NUM: u64 = 16;
