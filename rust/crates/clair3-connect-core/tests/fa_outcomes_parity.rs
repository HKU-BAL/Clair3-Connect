use clair3_connect_core::vcf::probabilities::{possible_outcomes_fa, Outcome};
use std::path::PathBuf;

const Y_COLS: usize = 90;
const OUT_COLS: usize = 26;

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
    let major = bytes[6];
    let header_len = if major == 1 {
        u16::from_le_bytes([bytes[8], bytes[9]]) as usize + 10
    } else {
        u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize + 12
    };
    let data = &bytes[header_len..];
    let total = data.len() / 4;
    let rows = total / expected_cols;
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
fn possible_outcomes_fa_matches_python_aggregated() {
    let dir = fixture_dir();
    let y_path = dir.join("Y_full_alignment_cpu.npy");
    let outcomes_path = dir.join("fa_outcomes.npy");
    let alt_path = dir.join("fa_haplotag_alt_info.tsv");
    if !y_path.exists() || !outcomes_path.exists() || !alt_path.exists() {
        eprintln!("Skipping — fixtures missing");
        return;
    }

    let (n_rows, y_flat) = read_npy_f32(&y_path, Y_COLS);
    let (n_o, expected) = read_npy_f32(&outcomes_path, OUT_COLS);
    assert_eq!(n_rows, n_o, "Y rows ({n_rows}) vs outcomes rows ({n_o})");
    let alt_text = std::fs::read_to_string(&alt_path).unwrap();
    let alt_lines: Vec<&str> = alt_text.lines().collect();
    assert_eq!(n_rows, alt_lines.len(), "alt_info rows differ");

    let mut max_abs = 0.0_f32;
    let mut row_with_max: Option<usize> = None;
    let mut col_with_max: Option<usize> = None;
    let mut shortcut_count = 0;
    let mut full_count = 0;

    for row in 0..n_rows {
        let gt21 = &y_flat[row * Y_COLS..row * Y_COLS + 21];
        let geno = &y_flat[row * Y_COLS + 21..row * Y_COLS + 24];
        let vlen1 = &y_flat[row * Y_COLS + 24..row * Y_COLS + 57];
        let vlen2 = &y_flat[row * Y_COLS + 57..row * Y_COLS + 90];

        let alt_line = alt_lines[row];
        let parts: Vec<&str> = alt_line.splitn(4, '-').collect();
        let ref_char = parts
            .get(2)
            .and_then(|s| s.bytes().next())
            .unwrap_or(b'N');

        let out = possible_outcomes_fa(gt21, geno, vlen1, vlen2, ref_char);

        let mut got = [0.0_f32; OUT_COLS];
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
                got[12] = f.homo_ins_probabilities.iter().sum();
                got[13] = f.hetero_insins_probabilities.iter().sum();

                for k in 0..4 {
                    got[14 + k] = f
                        .hetero_acgt_ins_probabilities
                        .iter()
                        .skip(k)
                        .step_by(4)
                        .sum();
                }
                got[18] = f.homo_del_probabilities.iter().sum();
                got[19] = f.hetero_deldel_probabilities.iter().sum();
                for k in 0..4 {
                    got[20 + k] = f
                        .hetero_acgt_del_probabilities
                        .iter()
                        .skip(k)
                        .step_by(4)
                        .sum();
                }
                got[24] = f.hetero_insdel_probabilities.iter().sum();
            }
        }

        for col in 0..OUT_COLS {
            let g = got[col];
            let e = expected[row * OUT_COLS + col];
            if !g.is_finite() && !e.is_finite() {
                continue;
            }
            let d = (g - e).abs();
            if d > max_abs {
                max_abs = d;
                row_with_max = Some(row);
                col_with_max = Some(col);
            }
        }
    }

    eprintln!(
        "rows={n_rows} (shortcut={shortcut_count} full={full_count}) \
         max_abs_diff={max_abs:e} at row={row_with_max:?} col={col_with_max:?}"
    );

    assert!(
        max_abs < 1e-5,
        "FA outcomes diverged: max_abs={max_abs:e} (col={col_with_max:?})"
    );
}
