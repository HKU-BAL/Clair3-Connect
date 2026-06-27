use clair3_connect_core::vcf::output_with::{output_with, OutputWithConfig, OutputWithInputs};
use std::path::PathBuf;

const Y_COLS: usize = 24;

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
    assert_eq!(&bytes[0..6], b"\x93NUMPY");
    let major = bytes[6];
    let header_len = if major == 1 {
        u16::from_le_bytes([bytes[8], bytes[9]]) as usize + 10
    } else {
        u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize + 12
    };
    let header = std::str::from_utf8(&bytes[10..header_len]).unwrap();
    assert!(header.contains("'<f4'") || header.contains("'float32'"));
    let data = &bytes[header_len..];
    let total = data.len() / 4;
    let rows = total / expected_cols;
    assert_eq!(rows * expected_cols * 4, data.len());
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
fn pileup_vcf_rows_match_python() {
    let dir = fixture_dir();
    let y_path = dir.join("Y_pileup_cpu.npy");
    let pos_alt_path = dir.join("pileup_windowed_pos_alt.tsv");
    let python_vcf = dir.join("pileup.vcf");

    if !y_path.exists() || !python_vcf.exists() {
        eprintln!("Skipping — Phase 0 fixture not present");
        return;
    }

    let (n_rows, y_flat) = read_npy_f32(&y_path, Y_COLS);
    let pos_alt = std::fs::read_to_string(&pos_alt_path).expect("read pos_alt");
    let pos_alt_lines: Vec<&str> = pos_alt.lines().collect();
    assert_eq!(n_rows, pos_alt_lines.len());

    let config = OutputWithConfig {
        is_show_reference: true,
        ..Default::default()
    };

    let mut got = String::new();
    for row in 0..n_rows {
        let (chr_pos_ref, payload) = pos_alt_lines[row].split_once('\t').unwrap();
        let chr_pos_seq = chr_pos_ref;
        let gt21 = &y_flat[row * Y_COLS..row * Y_COLS + 21];
        let geno = &y_flat[row * Y_COLS + 21..row * Y_COLS + 24];
        let res = output_with(
            OutputWithInputs {
                chr_pos_seq,
                alt_info: payload,
                gt21_probabilities: gt21,
                genotype_probabilities: geno,
                vlen_heads: None,
            },
            &config,
        )
        .expect("output_with");
        if let Some(line) = res {
            got.push_str(&line);
        }
    }

    let python_text = std::fs::read_to_string(&python_vcf).expect("read python vcf");
    let python_rows: Vec<&str> = python_text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect();
    let rust_rows: Vec<&str> = got.lines().collect();

    eprintln!("python_rows={} rust_rows={}", python_rows.len(), rust_rows.len());

    let mut mismatches: Vec<(usize, String, String)> = Vec::new();
    let max_check = python_rows.len().min(rust_rows.len());
    for i in 0..max_check {
        if python_rows[i] != rust_rows[i] {
            mismatches.push((i, python_rows[i].to_string(), rust_rows[i].to_string()));
            if mismatches.len() >= 5 {
                break;
            }
        }
    }
    if !mismatches.is_empty() {
        for (i, p, r) in &mismatches {
            eprintln!("row {i}:\n  python: {p}\n  rust:   {r}");
        }
    }
    assert_eq!(python_rows.len(), rust_rows.len(), "row count differs");
    assert!(
        mismatches.is_empty(),
        "{} row mismatches; first 5 above",
        mismatches.len()
    );
}
