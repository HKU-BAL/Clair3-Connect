use clair3_connect_core::vcf::probabilities::{possible_outcomes_pileup, Outcome};
use std::path::PathBuf;

const Y_COLS: usize = 24;
const OUTCOMES_COLS: usize = 26;

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("test_data/ont_chr20_demo/artifacts")
}

fn read_npy_f32(path: &std::path::Path, expected_cols: usize) -> (usize, Vec<f32>) {

    let bytes = std::fs::read(path).expect("read npy");

    assert_eq!(&bytes[0..6], b"\x93NUMPY", "not a npy file: {:?}", path);
    let major = bytes[6];
    let _minor = bytes[7];
    let header_len = if major == 1 {
        u16::from_le_bytes([bytes[8], bytes[9]]) as usize + 10
    } else {
        u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize + 12
    };
    let header = std::str::from_utf8(&bytes[10..header_len]).unwrap();
    assert!(
        header.contains("'<f4'") || header.contains("'float32'"),
        "expected float32 dtype: header={header}"
    );
    let data = &bytes[header_len..];
    let total = data.len() / 4;
    let rows = total / expected_cols;
    assert_eq!(rows * expected_cols * 4, data.len(), "shape mismatch");
    let mut out = Vec::with_capacity(total);
    for i in 0..total {
        let off = i * 4;
        out.push(f32::from_le_bytes([
            data[off],
            data[off + 1],
            data[off + 2],
            data[off + 3],
        ]));
    }
    (rows, out)
}

#[test]
fn pileup_outcomes_bit_equal_to_python() {
    let dir = fixture_dir();
    let y_path = dir.join("Y_pileup_cpu.npy");
    let outcomes_path = dir.join("pileup_outcomes.npy");
    let pos_alt_path = dir.join("pileup_windowed_pos_alt.tsv");

    if !y_path.exists() || !outcomes_path.exists() {
        eprintln!("Skipping — Phase 0 fixture not present at {:?}", dir);
        return;
    }

    let (n_rows, y_flat) = read_npy_f32(&y_path, Y_COLS);
    let (n_rows_o, expected_flat) = read_npy_f32(&outcomes_path, OUTCOMES_COLS);
    assert_eq!(n_rows, n_rows_o, "row count mismatch Y vs outcomes");

    let pos_alt_text = std::fs::read_to_string(&pos_alt_path).expect("read pos_alt");
    let pos_alt_lines: Vec<&str> = pos_alt_text.lines().collect();
    assert_eq!(
        pos_alt_lines.len(),
        n_rows,
        "pos_alt rows differ from Y rows"
    );

    let mut max_abs_diff = 0.0_f32;
    let mut row_with_max: Option<usize> = None;
    let mut col_with_max: Option<usize> = None;
    let mut shortcut_count = 0;
    let mut full_count = 0;

    for row in 0..n_rows {
        let gt21 = &y_flat[row * Y_COLS..row * Y_COLS + 21];
        let geno = &y_flat[row * Y_COLS + 21..row * Y_COLS + 24];
        let line = pos_alt_lines[row];
        let chr_pos_ref = line.split('\t').next().unwrap();
        let ref_base = chr_pos_ref.rsplit(':').next().unwrap().as_bytes()[0];

        let out = possible_outcomes_pileup(gt21, geno, ref_base);

        let mut got = [0.0_f32; OUTCOMES_COLS];
        match out {
            Outcome::HomoRefShortcut(p) => {
                shortcut_count += 1;
                got[0] = 1.0;
                got[1] = p;
            }
            Outcome::Full(f) => {
                full_count += 1;
                got[0] = 0.0;
                got[1] = f.homo_ref;
                got[2..6].copy_from_slice(&f.homo_snp);
                got[6..12].copy_from_slice(&f.hetero_snp);
                got[12] = f.homo_ins_probabilities[0];
                got[13] = f.hetero_insins_probabilities[0];
                got[14..18].copy_from_slice(&f.hetero_acgt_ins_probabilities);
                got[18] = f.homo_del_probabilities[0];
                got[19] = f.hetero_deldel_probabilities[0];
                got[20..24].copy_from_slice(&f.hetero_acgt_del_probabilities);
                got[24] = f.hetero_insdel_probabilities[0];
            }
        }

        for col in 0..OUTCOMES_COLS {
            let exp = expected_flat[row * OUTCOMES_COLS + col];
            let g = got[col];
            assert!(
                g.is_finite() && exp.is_finite() || (g.is_nan() == exp.is_nan()),
                "non-finite mismatch row={row} col={col} got={g} exp={exp}"
            );
            let d = (g - exp).abs();
            if d > max_abs_diff {
                max_abs_diff = d;
                row_with_max = Some(row);
                col_with_max = Some(col);
            }
        }
    }

    eprintln!(
        "rows={n_rows} (shortcut={shortcut_count} full={full_count}) \
         max_abs_diff={max_abs_diff:e} at row={row_with_max:?} col={col_with_max:?}",
    );
    assert_eq!(
        max_abs_diff, 0.0,
        "Rust outcomes diverge from Python (max_abs={max_abs_diff:e} at \
         row={row_with_max:?} col={col_with_max:?})"
    );
}
