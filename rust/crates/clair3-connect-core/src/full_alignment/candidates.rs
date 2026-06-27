use super::constants::FLANKING_BASE_NUM;
use std::collections::HashMap;

#[derive(Debug)]
pub struct FaCandidateContext {

    pub candidates: Vec<u64>,

    pub candidate_index: HashMap<u64, usize>,

    pub flanking_index: HashMap<u64, usize>,

    pub flanking_positions: Vec<u64>,
}

impl FaCandidateContext {

    pub fn build(candidates: &[u64]) -> Self {
        debug_assert!(
            candidates.windows(2).all(|w| w[0] < w[1]),
            "candidates must be sorted and unique"
        );

        let candidate_index: HashMap<u64, usize> = candidates
            .iter()
            .enumerate()
            .map(|(i, p)| (*p, i))
            .collect();

        let mut flanking_index: HashMap<u64, usize> = HashMap::new();
        let mut flanking_positions: Vec<u64> = Vec::new();
        for c in candidates {
            for j in 0..(2 * FLANKING_BASE_NUM + 1) {

                let signed = *c as i64 - FLANKING_BASE_NUM as i64 + j as i64;
                let key = if signed < 0 { 0 } else { signed as u64 };
                flanking_index.entry(key).or_insert_with(|| {
                    let idx = flanking_positions.len();
                    flanking_positions.push(key);
                    idx
                });
            }
        }

        Self {
            candidates: candidates.to_vec(),
            candidate_index,
            flanking_index,
            flanking_positions,
        }
    }

    pub fn candidate_num(&self) -> usize {
        self.candidates.len()
    }

    pub fn flanking_num(&self) -> usize {
        self.flanking_positions.len()
    }

    pub fn flanking_idx(&self, pos: u64) -> Option<usize> {
        self.flanking_index.get(&pos).copied()
    }

    pub fn candidate_idx(&self, pos: u64) -> Option<usize> {
        self.candidate_index.get(&pos).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_candidate_produces_33_flanking_positions() {
        let ctx = FaCandidateContext::build(&[100]);
        assert_eq!(ctx.candidate_num(), 1);
        assert_eq!(ctx.flanking_num(), 33);

        assert_eq!(ctx.flanking_idx(84), Some(0));
        assert_eq!(ctx.flanking_idx(100), Some(16));
        assert_eq!(ctx.flanking_idx(116), Some(32));
        assert_eq!(ctx.flanking_idx(83), None);
        assert_eq!(ctx.flanking_idx(117), None);
    }

    #[test]
    fn overlapping_candidates_dedup_flanking() {

        let ctx = FaCandidateContext::build(&[100, 110]);
        assert_eq!(ctx.candidate_num(), 2);
        assert_eq!(ctx.flanking_num(), 43);
    }

    #[test]
    fn candidate_near_zero_clamps_to_zero() {

        let ctx = FaCandidateContext::build(&[10]);
        assert_eq!(ctx.flanking_idx(0), Some(0));

        assert_eq!(ctx.flanking_num(), 27);
        assert_eq!(ctx.flanking_idx(26), Some(26));
        assert_eq!(ctx.flanking_idx(27), None);
    }

    #[test]
    fn lookups_correct() {
        let ctx = FaCandidateContext::build(&[100, 200]);
        assert_eq!(ctx.candidate_idx(100), Some(0));
        assert_eq!(ctx.candidate_idx(200), Some(1));
        assert_eq!(ctx.candidate_idx(150), None);
        assert_eq!(ctx.flanking_idx(100), Some(16));
    }
}
