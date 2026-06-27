use std::fs::File;
use std::path::PathBuf;

use ndarray::{Array1, Array2};
use ndarray_npy::ReadNpyExt;

use clair3_connect_core::hts::region::Region;
use clair3_connect_core::pileup::kernel::{calculate_pileup, PileupParams};

const _FIXTURE_REGION_USER: &str = "chr20:100000-300000";

const FIXTURE_REGION_PRODUCER_EFFECTIVE: &str = "chr20:99966-300033";

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

    let region = Region::parse(FIXTURE_REGION_PRODUCER_EFFECTIVE).expect("parse region");
    let params = PileupParams::default();
    let res = calculate_pileup(&region, &bam, &ref_, &params).expect("calculate_pileup");
    Some(res)
}

#[test]
fn m1_column_count() {
    let Some(dir) = golden_dir() else {
        eprintln!("[skip] CLAIR3_CONNECT_GOLDEN_DIR not set");
        return;
    };
    let major_path = dir.join("pileup_raw_chunk_0000_major.npy");
    if !major_path.exists() {
        eprintln!("[skip] {} not present", major_path.display());
        return;
    }
    let golden_major =
        Array1::<i64>::read_npy(File::open(&major_path).unwrap()).expect("read major npy");
    let expected_n_cols = golden_major.len();

    let Some(res) = run_fixture_pileup() else {
        return;
    };

    assert_eq!(
        res.chunks.len(),
        1,
        "expected one chunk on no-gap fixture, got {}",
        res.chunks.len()
    );
    let chunk = &res.chunks[0];
    eprintln!(
        "[m1] rust n_cols = {}, golden n_cols = {}, chunks = {}",
        chunk.n_cols,
        expected_n_cols,
        res.chunks.len()
    );
    assert_eq!(
        chunk.n_cols, expected_n_cols,
        "column count mismatch: rust {} vs golden {}",
        chunk.n_cols, expected_n_cols
    );
}

#[test]
fn m2_major_minor_match() {
    let Some(dir) = golden_dir() else {
        eprintln!("[skip] CLAIR3_CONNECT_GOLDEN_DIR not set");
        return;
    };
    let major_path = dir.join("pileup_raw_chunk_0000_major.npy");
    let minor_path = dir.join("pileup_raw_chunk_0000_minor.npy");
    if !major_path.exists() || !minor_path.exists() {
        eprintln!("[skip] missing golden npy");
        return;
    }
    let golden_major = Array1::<i64>::read_npy(File::open(&major_path).unwrap()).unwrap();
    let golden_minor = Array1::<i64>::read_npy(File::open(&minor_path).unwrap()).unwrap();

    let Some(res) = run_fixture_pileup() else {
        return;
    };
    let chunk = &res.chunks[0];
    assert_eq!(chunk.n_cols, golden_major.len(), "n_cols mismatch");

    for i in 0..chunk.n_cols {
        let r = chunk.major[i];
        let g = golden_major[i];
        if r != g {

            let ctx_start = i.saturating_sub(3);
            let ctx_end = (i + 4).min(chunk.n_cols);
            for k in ctx_start..ctx_end {
                eprintln!(
                    "  [{}] rust={} golden={}",
                    k, chunk.major[k], golden_major[k]
                );
            }
            panic!("major[{}] mismatch: rust {} vs golden {}", i, r, g);
        }
    }
    for i in 0..chunk.n_cols {
        assert_eq!(
            chunk.minor[i], golden_minor[i],
            "minor[{}] mismatch: rust {} vs golden {}",
            i, chunk.minor[i], golden_minor[i]
        );
    }
    eprintln!("[m2] major/minor match across {} cols", chunk.n_cols);
}

#[test]
fn m3_base_counts_match() {
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

    let cols = [0_usize, 1, 2, 3, 9, 10, 11, 12, 8, 17];
    let mut first_mismatch: Option<(usize, usize, i64, i64)> = None;
    let mut mismatch_count = 0usize;
    for row in 0..chunk.n_cols {
        for &c in &cols {
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
        eprintln!("[m3] first mismatch at row={row} pos={pos} col={c}");
        eprintln!("     rust   row: {rust_row:?}");
        eprintln!("     golden row: {golden_row:?}");
        eprintln!("     rust = {r}, golden = {g}");
        eprintln!("[m3] total mismatches: {mismatch_count}");
    }
    assert_eq!(mismatch_count, 0, "M3 column mismatches: {mismatch_count}");
    eprintln!("[m3] base counts match across {} cols", chunk.n_cols);
}

#[test]
fn m4_ins_del_buffers_match() {
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

    let cols = [4_usize, 5, 6, 7, 13, 14, 15, 16];
    let mut first_mismatch: Option<(usize, usize, i64, i64)> = None;
    let mut mismatch_count = 0usize;
    for row in 0..chunk.n_cols {
        for &c in &cols {
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
        eprintln!("[m4] first mismatch at row={row} pos={pos} col={c}");
        eprintln!("     rust   row: {rust_row:?}");
        eprintln!("     golden row: {golden_row:?}");
        eprintln!("     rust cell = {r}, golden cell = {g}");
        eprintln!("[m4] total mismatches: {mismatch_count}");
    }
    assert_eq!(
        mismatch_count, 0,
        "M4 ins/del summary mismatches: {mismatch_count}"
    );
    eprintln!(
        "[m4] ins/del summary cols match across {} cols × {} cells",
        chunk.n_cols,
        cols.len()
    );
}
