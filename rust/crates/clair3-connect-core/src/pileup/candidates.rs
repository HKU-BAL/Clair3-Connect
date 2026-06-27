use super::kernel::{AltInfo, PileupChunk};
use super::tables::{FEAT_LEN, FLANKING_BASE_NUM};

#[derive(Debug, Clone)]
pub struct WindowedCandidate {

    pub tensor: Vec<i32>,

    pub pos_info: String,

    pub alt_info: String,
}

const NO_OF_POSITIONS: usize = (2 * FLANKING_BASE_NUM + 1) as usize;

pub fn window_candidates(chunks: &[PileupChunk], candidates: &[AltInfo]) -> Vec<WindowedCandidate> {
    let mut out: Vec<WindowedCandidate> = Vec::with_capacity(candidates.len());
    for cand in candidates {

        let pos_0based: i64 = cand.pos_1based - 1;
        let start = pos_0based - FLANKING_BASE_NUM as i64;

        let end_for_python_check = pos_0based + FLANKING_BASE_NUM as i64 + 2;

        let mut emitted = false;
        for chunk in chunks {
            if chunk.n_cols == 0 {
                continue;
            }
            let first = chunk.major[0];
            let last = chunk.major[chunk.n_cols - 1];
            if start >= first && end_for_python_check <= last {

                let offset = (start - first) as usize;

                if offset + NO_OF_POSITIONS > chunk.n_cols {
                    break;
                }

                let slice_start = offset * FEAT_LEN;
                let slice_end = slice_start + NO_OF_POSITIONS * FEAT_LEN;
                let raw = &chunk.matrix[slice_start..slice_end];

                let mut any_zero_col = false;
                for r in 0..NO_OF_POSITIONS {
                    let row_start = r * FEAT_LEN;
                    let row_end = row_start + FEAT_LEN;
                    if raw[row_start..row_end].iter().all(|&v| v == 0) {
                        any_zero_col = true;
                        break;
                    }
                }
                if any_zero_col {
                    break;
                }

                let tensor: Vec<i32> = raw.iter().map(|&v| v as i32).collect();
                out.push(WindowedCandidate {
                    tensor,
                    pos_info: cand.pos_info.clone(),
                    alt_info: cand.alt_info.clone(),
                });
                emitted = true;
                break;
            }
        }
        let _ = emitted;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::kernel::PileupChunk;
    use crate::pileup::tables::FEAT_LEN;

    fn make_chunk(start: i64, n: usize) -> PileupChunk {
        let mut major: Vec<i64> = Vec::with_capacity(n);
        for i in 0..n {
            major.push(start + i as i64);
        }
        let matrix: Vec<i64> = (0..n * FEAT_LEN).map(|i| (i + 1) as i64).collect();
        let minor = vec![0i64; n];
        PileupChunk {
            matrix,
            major,
            minor,
            n_cols: n,
        }
    }

    fn alt_info_at(pos_1based: i64) -> AltInfo {
        AltInfo {
            pos_1based,
            pos_info: format!("chr1:{}:A", pos_1based),
            alt_info: "10-XC 5".to_string(),
        }
    }

    #[test]
    fn candidate_in_middle_emits_full_window() {

        let chunk = make_chunk(100, 100);
        let cand = alt_info_at(150);
        let out = window_candidates(&[chunk], &[cand]);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].tensor.len(), NO_OF_POSITIONS * FEAT_LEN);
        assert_eq!(out[0].pos_info, "chr1:150:A");

        let expected_first = (33 * FEAT_LEN + 1) as i32;
        assert_eq!(out[0].tensor[0], expected_first);
    }

    #[test]
    fn candidate_too_close_to_chunk_start_is_dropped() {
        let chunk = make_chunk(100, 100);

        let cand = alt_info_at(110);
        let out = window_candidates(&[chunk], &[cand]);
        assert!(out.is_empty());
    }

    #[test]
    fn candidate_too_close_to_chunk_end_is_dropped() {
        let chunk = make_chunk(100, 50);

        let cand = alt_info_at(145);
        let out = window_candidates(&[chunk], &[cand]);
        assert!(out.is_empty());
    }

    #[test]
    fn right_boundary_matches_python_strictness() {
        let chunk = make_chunk(100, 100);
        let chunks = std::slice::from_ref(&chunk);

        let out = window_candidates(chunks, &[alt_info_at(182)]);
        assert_eq!(out.len(), 1, "pos_1based=182 should be the last accepted");

        let out = window_candidates(chunks, &[alt_info_at(183)]);
        assert!(
            out.is_empty(),
            "pos_1based=183 should be rejected (matches Python)"
        );
    }

    #[test]
    fn all_zero_column_drops_candidate() {
        let mut chunk = make_chunk(100, 100);

        let zero_offset = 33;
        for c in 0..FEAT_LEN {
            chunk.matrix[zero_offset * FEAT_LEN + c] = 0;
        }
        let cand = alt_info_at(150);
        let out = window_candidates(&[chunk], &[cand]);
        assert!(
            out.is_empty(),
            "expected candidate dropped due to all-zero col"
        );
    }

    #[test]
    fn multiple_candidates_in_one_chunk() {
        let chunk = make_chunk(100, 200);
        let cands = vec![alt_info_at(150), alt_info_at(160), alt_info_at(170)];
        let out = window_candidates(&[chunk], &cands);
        assert_eq!(out.len(), 3);
    }
}
