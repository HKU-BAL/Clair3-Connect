use super::kernel::PileupChunk;
use super::tables::FEAT_LEN;

pub fn split_on_gaps(chunk: PileupChunk) -> Vec<PileupChunk> {
    let n = chunk.major.len();
    debug_assert_eq!(chunk.n_cols, n);
    debug_assert_eq!(chunk.minor.len(), n);
    debug_assert_eq!(chunk.matrix.len(), n * FEAT_LEN);

    debug_assert!(
        chunk.major.windows(2).all(|w| w[0] < w[1]),
        "split_on_gaps requires strictly increasing `major`"
    );

    if n == 0 {
        return Vec::new();
    }

    let has_gap = (1..n).any(|i| chunk.major[i] - chunk.major[i - 1] > 1);
    if !has_gap {
        return vec![chunk];
    }

    let mut splits: Vec<usize> = Vec::with_capacity(8);
    splits.push(0);
    for i in 1..n {
        if chunk.major[i] - chunk.major[i - 1] > 1 {
            splits.push(i);
        }
    }
    splits.push(n);

    let mut out: Vec<PileupChunk> = Vec::with_capacity(splits.len() - 1);
    for w in splits.windows(2) {
        let (start, end) = (w[0], w[1]);
        if start == end {
            continue;
        }
        let matrix = chunk.matrix[start * FEAT_LEN..end * FEAT_LEN].to_vec();
        let major = chunk.major[start..end].to_vec();
        let minor = chunk.minor[start..end].to_vec();
        out.push(PileupChunk {
            matrix,
            major,
            minor,
            n_cols: end - start,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pileup::kernel::PileupChunk;
    use crate::pileup::tables::FEAT_LEN;

    fn make_chunk(positions: &[i64]) -> PileupChunk {
        let n = positions.len();
        let matrix = (0..n * FEAT_LEN).map(|i| i as i64).collect();
        let major = positions.to_vec();
        let minor = vec![0i64; n];
        PileupChunk {
            matrix,
            major,
            minor,
            n_cols: n,
        }
    }

    #[test]
    fn empty_in_empty_out() {
        let c = make_chunk(&[]);
        assert!(split_on_gaps(c).is_empty());
    }

    #[test]
    fn single_column() {
        let c = make_chunk(&[42]);
        let out = split_on_gaps(c);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].major, vec![42]);
        assert_eq!(out[0].matrix.len(), FEAT_LEN);
    }

    #[test]
    fn no_gap_returns_single_chunk() {
        let c = make_chunk(&[100, 101, 102, 103, 104]);
        let out = split_on_gaps(c);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].n_cols, 5);
        assert_eq!(out[0].major, vec![100, 101, 102, 103, 104]);
        assert_eq!(out[0].matrix.len(), 5 * FEAT_LEN);
    }

    #[test]
    fn single_gap_splits_in_two() {

        let c = make_chunk(&[100, 101, 102, 110, 111]);
        let out = split_on_gaps(c);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].major, vec![100, 101, 102]);
        assert_eq!(out[1].major, vec![110, 111]);

        assert_eq!(out[0].matrix.len(), 3 * FEAT_LEN);
        assert_eq!(out[1].matrix.len(), 2 * FEAT_LEN);

        assert_eq!(out[1].matrix[0], (3 * FEAT_LEN) as i64);
    }

    #[test]
    fn diff_of_one_is_not_a_gap() {

        let c = make_chunk(&[100, 101, 102]);
        let out = split_on_gaps(c);
        assert_eq!(out.len(), 1);
    }

    #[test]
    fn diff_of_two_is_a_gap() {

        let c = make_chunk(&[100, 102]);
        let out = split_on_gaps(c);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].major, vec![100]);
        assert_eq!(out[1].major, vec![102]);
    }

    #[test]
    fn multi_gap_splits() {
        let c = make_chunk(&[10, 11, 20, 21, 22, 50]);
        let out = split_on_gaps(c);
        assert_eq!(out.len(), 3);
        assert_eq!(out[0].major, vec![10, 11]);
        assert_eq!(out[1].major, vec![20, 21, 22]);
        assert_eq!(out[2].major, vec![50]);
    }

    #[test]
    fn matrix_n_cols_minor_lengths_preserved() {
        let c = make_chunk(&[100, 200, 201, 300, 301, 302]);
        let out = split_on_gaps(c);
        assert_eq!(out.len(), 3);
        for chunk in &out {
            assert_eq!(chunk.n_cols, chunk.major.len());
            assert_eq!(chunk.minor.len(), chunk.major.len());
            assert_eq!(chunk.matrix.len(), chunk.major.len() * FEAT_LEN);
        }
    }
}
