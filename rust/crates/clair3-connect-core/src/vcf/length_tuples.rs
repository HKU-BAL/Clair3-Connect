use super::labels::variant_length;

pub fn homo_ins_tuples(
    vlen1: &[f32],
    vlen2: &[f32],
    extra_probability: f32,
) -> Vec<(i32, f32)> {
    let offset = variant_length::INDEX_OFFSET;
    let max = variant_length::MAX;
    (1..=max)
        .map(|i| {
            let idx = (i as usize) + offset;
            (i, vlen1[idx] * vlen2[idx] * extra_probability)
        })
        .collect()
}

pub fn hetero_ins_tuples(vlen1: &[f32], vlen2: &[f32]) -> Vec<(i32, f32)> {
    let offset = variant_length::INDEX_OFFSET;
    let max = variant_length::MAX;
    (1..=max)
        .map(|i| {
            let idx = (i as usize) + offset;
            (i, vlen1[offset] * vlen2[idx])
        })
        .collect()
}

pub fn hetero_insins_tuples(
    vlen1: &[f32],
    vlen2: &[f32],
    extra_probability: f32,
) -> Vec<((i32, i32), f32)> {
    let offset = variant_length::INDEX_OFFSET;
    let max = variant_length::MAX;
    let mut out = Vec::new();
    for i in 1..=max {
        for j in i..=max {
            let p = vlen1[(i as usize) + offset]
                * vlen2[(j as usize) + offset]
                * extra_probability;
            out.push(((i, j), p));
        }
    }
    out
}

pub fn homo_del_tuples(
    vlen1: &[f32],
    vlen2: &[f32],
    extra_probability: f32,
) -> Vec<(i32, f32)> {
    let offset = variant_length::INDEX_OFFSET as i32;
    let max = variant_length::MAX;
    (1..=max)
        .map(|i| {
            let idx = (-i + offset) as usize;
            (i, vlen1[idx] * vlen2[idx] * extra_probability)
        })
        .collect()
}

pub fn hetero_del_tuples(vlen1: &[f32], vlen2: &[f32]) -> Vec<(i32, f32)> {
    let offset = variant_length::INDEX_OFFSET;
    let max = variant_length::MAX;
    (1..=max)
        .map(|i| {
            let neg_idx = (-i + offset as i32) as usize;
            (i, vlen1[neg_idx] * vlen2[offset])
        })
        .collect()
}

pub fn hetero_deldel_tuples(
    vlen1: &[f32],
    vlen2: &[f32],
    extra_probability: f32,
) -> Vec<((i32, i32), f32)> {
    let offset = variant_length::INDEX_OFFSET as i32;
    let max = variant_length::MAX;
    let mut out = Vec::new();
    for i in 1..=max {
        for j in 1..=max {
            if i == j && i as usize != variant_length::INDEX_OFFSET
                && j as usize != variant_length::INDEX_OFFSET
            {
                continue;
            }
            let tup = if i < j { (i, j) } else { (j, i) };
            let p = vlen1[(-i + offset) as usize]
                * vlen2[(-j + offset) as usize]
                * extra_probability;
            out.push((tup, p));
        }
    }
    out
}

pub fn hetero_insdel_tuples(
    vlen1: &[f32],
    vlen2: &[f32],
    extra_probability: f32,
) -> Vec<((i32, i32), f32)> {
    let offset = variant_length::INDEX_OFFSET as i32;
    let max = variant_length::MAX;
    let mut out = Vec::new();
    for i in 1..=max {
        for j in 1..=max {
            let p = vlen1[(-i + offset) as usize]
                * vlen2[(j as usize) + offset as usize]
                * extra_probability;
            out.push(((i, j), p));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ones_vlen() -> Vec<f32> {
        vec![1.0_f32; variant_length::OUTPUT_LABEL_COUNT]
    }

    #[test]
    fn homo_ins_returns_max_entries() {
        let v = homo_ins_tuples(&ones_vlen(), &ones_vlen(), 0.5);
        assert_eq!(v.len(), variant_length::MAX as usize);
        for (i, p) in &v {
            assert!(*i >= 1);
            assert!((*p - 0.5).abs() < 1e-6);
        }
    }

    #[test]
    fn hetero_insins_yields_upper_triangle() {

        let v = hetero_insins_tuples(&ones_vlen(), &ones_vlen(), 0.1);
        assert_eq!(v.len(), 16 * 17 / 2);
        for ((i, j), _) in &v {
            assert!(i <= j);
            assert!(*i >= 1 && *j <= 16);
        }
    }

    #[test]
    fn homo_del_uses_negative_offset() {
        let mut v = ones_vlen();

        v[15] = 9.0;
        let r = homo_del_tuples(&v, &v, 1.0);

        assert_eq!(r[0].0, 1);
        assert!((r[0].1 - 81.0).abs() < 1e-3);
    }

    #[test]
    fn hetero_insdel_combines_neg_and_pos_indices() {
        let v = ones_vlen();
        let r = hetero_insdel_tuples(&v, &v, 1.0);
        assert_eq!(r.len(), 16 * 16);
    }

    #[test]
    fn hetero_deldel_skips_same_length() {

        let r = hetero_deldel_tuples(&ones_vlen(), &ones_vlen(), 1.0);
        assert_eq!(r.len(), 16 * 16 - 15);
    }

    #[test]
    fn homo_ins_index_offset_arithmetic() {

        let mut v = vec![1.0_f32; 33];
        v[17] = 5.0;
        let r = homo_ins_tuples(&v, &v, 1.0);

        assert_eq!(r[0].0, 1);
        assert!((r[0].1 - 25.0).abs() < 1e-3);
    }
}
