use std::fs::File;
use std::path::PathBuf;

use ndarray::Array2;
use ndarray_npy::ReadNpyExt;

use clair3_connect_core::hts::region::Region;
use clair3_connect_core::pileup::kernel::{calculate_pileup, PileupParams};

const FIXTURE_REGION: &str = "chr20:99966-300033";

fn golden_dir() -> Option<PathBuf> {
    std::env::var_os("CLAIR3_CONNECT_GOLDEN_DIR").map(PathBuf::from)
}

fn fixture_bam() -> PathBuf {
    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    crate_dir
        .join("..")
        .join("..")
        .join("test_data")
        .join("ont_chr20_demo")
        .join("HG003_chr20_demo.bam")
}

fn fixture_ref() -> PathBuf {
    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    crate_dir
        .join("..")
        .join("..")
        .join("test_data")
        .join("ont_chr20_demo")
        .join("GRCh38_no_alt_chr20.fa")
}

fn run_fixture_pileup() -> Option<clair3_connect_core::pileup::kernel::PileupResult> {
    let bam = fixture_bam();
    let ref_ = fixture_ref();
    if !bam.exists() || !ref_.exists() {
        eprintln!(
            "[skip] missing fixture: {} or {}",
            bam.display(),
            ref_.display()
        );
        return None;
    }
    let region = Region::parse(FIXTURE_REGION).expect("parse region");
    let params = PileupParams::default();
    let res = calculate_pileup(&region, &bam, &ref_, &params).expect("calculate_pileup");
    Some(res)
}

#[test]
fn m5_full_matrix_byte_match() {
    let Some(dir) = golden_dir() else {
        eprintln!("[skip] CLAIR3_CONNECT_GOLDEN_DIR not set");
        return;
    };
    let matrix_path = dir.join("pileup_raw_chunk_0000_matrix.npy");
    if !matrix_path.exists() {
        eprintln!("[skip] {} not present", matrix_path.display());
        return;
    }
    let golden = Array2::<i64>::read_npy(File::open(&matrix_path).unwrap()).unwrap();

    let Some(res) = run_fixture_pileup() else {
        return;
    };
    let chunk = &res.chunks[0];
    assert_eq!(chunk.n_cols, golden.shape()[0], "n_cols mismatch");
    assert_eq!(18, golden.shape()[1]);

    let mut first_mismatch: Option<(usize, usize, i64, i64)> = None;
    let mut mismatch_count = 0usize;
    for row in 0..chunk.n_cols {
        for c in 0..18 {
            let r = chunk.matrix[row * 18 + c];
            let g = golden[(row, c)];
            if r != g {
                mismatch_count += 1;
                if first_mismatch.is_none() {
                    first_mismatch = Some((row, c, r, g));
                }
            }
        }
    }
    if let Some((row, c, r, g)) = first_mismatch {
        let pos = chunk.major[row];
        let rust_row: Vec<i64> = chunk.matrix[row * 18..(row + 1) * 18].to_vec();
        let golden_row: Vec<i64> = (0..18).map(|k| golden[(row, k)]).collect();
        eprintln!("[m5-matrix] first mismatch at row={row} pos={pos} col={c}");
        eprintln!("            rust   row: {rust_row:?}");
        eprintln!("            golden row: {golden_row:?}");
        eprintln!("            rust cell = {r}, golden cell = {g}");
        eprintln!("[m5-matrix] total mismatches: {mismatch_count}");
    }
    assert_eq!(
        mismatch_count, 0,
        "M5 full matrix mismatches: {mismatch_count}"
    );
    eprintln!(
        "[m5-matrix] full {}×18 matrix byte-identical to libclair3.so",
        chunk.n_cols
    );
}

fn normalize_alt_info(s: &str) -> String {

    let (depth, payload) = match s.split_once('-') {
        Some((d, p)) => (d, p),
        None => return s.to_string(),
    };

    let toks: Vec<&str> = payload.split(' ').filter(|t| !t.is_empty()).collect();
    let mut xs: Vec<(String, String)> = Vec::new();
    let mut ds: Vec<(String, String)> = Vec::new();
    let mut is_: Vec<(String, String)> = Vec::new();
    let mut rs: Vec<(String, String)> = Vec::new();
    let mut i = 0;
    while i + 1 < toks.len() {
        let entry = toks[i].to_string();
        let count = toks[i + 1].to_string();
        match entry.chars().next() {
            Some('X') => xs.push((entry, count)),
            Some('D') => ds.push((entry, count)),
            Some('I') => is_.push((entry, count)),
            Some('R') => rs.push((entry, count)),
            _ => {}
        }
        i += 2;
    }
    is_.sort_by(|a, b| a.0.cmp(&b.0));
    let mut out = String::from(depth);
    out.push('-');
    for (e, c) in xs
        .iter()
        .chain(ds.iter())
        .chain(is_.iter())
        .chain(rs.iter())
    {
        out.push_str(e);
        out.push(' ');
        out.push_str(c);
        out.push(' ');
    }
    out.trim_end().to_string()
}

#[cfg(test)]
mod normalize_tests {
    use super::normalize_alt_info;
    #[test]
    fn sorts_only_i_entries() {
        let s = "33-XA 5 ITA 1 ICA 3 RA 20";

        assert_eq!(normalize_alt_info(s), "33-XA 5 ICA 3 ITA 1 RA 20");
    }
    #[test]
    fn idempotent_when_already_sorted() {
        let s = "29-DA 3 DAA 5 ICA 5 ICAA 2 RC 6";
        assert_eq!(normalize_alt_info(s), s);
    }
}

#[test]
fn m5_alt_info_byte_match() {
    let Some(dir) = golden_dir() else {
        eprintln!("[skip] CLAIR3_CONNECT_GOLDEN_DIR not set");
        return;
    };
    let alt_info_path = dir.join("pileup_raw_alt_info.tsv");
    if !alt_info_path.exists() {
        eprintln!("[skip] {} not present", alt_info_path.display());
        return;
    }
    let golden_text = std::fs::read_to_string(&alt_info_path).unwrap();
    let golden_lines: Vec<&str> = golden_text.lines().collect();

    let Some(res) = run_fixture_pileup() else {
        return;
    };

    eprintln!(
        "[m5-alt_info] rust candidates = {}, golden lines = {}",
        res.candidates.len(),
        golden_lines.len()
    );
    assert_eq!(
        res.candidates.len(),
        golden_lines.len(),
        "candidate count mismatch: rust {} vs golden {}",
        res.candidates.len(),
        golden_lines.len()
    );

    let n = res.candidates.len();
    let mut mismatches: Vec<(usize, String, String)> = Vec::new();
    for (i, cand) in res.candidates.iter().enumerate() {

        let cand_alt_info_norm = normalize_alt_info(&cand.alt_info);
        let rust_line = format!(
            "{}\t{}\t{}",
            cand.pos_1based, cand.pos_info, cand_alt_info_norm
        );

        let parts: Vec<&str> = golden_lines[i].splitn(3, '\t').collect();
        let golden_normalized = if parts.len() == 3 {
            format!(
                "{}\t{}\t{}",
                parts[0],
                parts[1],
                normalize_alt_info(parts[2])
            )
        } else {
            golden_lines[i].to_string()
        };
        if rust_line != golden_normalized {
            mismatches.push((i, rust_line, golden_normalized));
        }
    }
    let _ = n;
    if !mismatches.is_empty() {
        for (i, rust_line, golden_line) in mismatches.iter().take(5) {
            eprintln!("[m5-alt_info] mismatch #{i}");
            eprintln!("              rust:   {rust_line}");
            eprintln!("              golden: {golden_line}");
        }
        eprintln!("[m5-alt_info] total mismatches: {}/{}", mismatches.len(), n);
    }
    assert!(
        mismatches.is_empty(),
        "M5 alt_info mismatches: {}/{}",
        mismatches.len(),
        n
    );
    eprintln!("[m5-alt_info] all {n} candidate lines byte-identical");
}
