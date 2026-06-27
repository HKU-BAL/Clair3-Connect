use std::fs::File;
use std::path::PathBuf;

use ndarray::Array3;
use ndarray_npy::ReadNpyExt;

use clair3_connect_core::hts::region::Region;
use clair3_connect_core::pileup::candidates::window_candidates;
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

#[test]
fn m7_windowed_tensor_byte_match() {
    let Some(dir) = golden_dir() else {
        eprintln!("[skip] CLAIR3_CONNECT_GOLDEN_DIR not set");
        return;
    };
    let win_path = dir.join("pileup_windowed_X_int32.npy");
    if !win_path.exists() {
        eprintln!("[skip] {} not present", win_path.display());
        return;
    }
    let bam = fixture_bam();
    let ref_ = fixture_ref();
    if !bam.exists() || !ref_.exists() {
        eprintln!("[skip] missing fixture BAM/ref");
        return;
    }

    let golden = Array3::<i32>::read_npy(File::open(&win_path).unwrap()).unwrap();
    let (n_golden, h, w) = (golden.shape()[0], golden.shape()[1], golden.shape()[2]);
    assert_eq!((h, w), (33, 18), "unexpected windowed shape");

    let region = Region::parse(FIXTURE_REGION).expect("parse region");
    let params = PileupParams::default();
    let res = calculate_pileup(&region, &bam, &ref_, &params).expect("calculate_pileup");
    let windowed = window_candidates(&res.chunks, &res.candidates);

    eprintln!(
        "[m7] rust windowed = {}, golden N = {}, raw candidates = {}",
        windowed.len(),
        n_golden,
        res.candidates.len()
    );
    assert_eq!(
        windowed.len(),
        n_golden,
        "windowed-candidate count mismatch: rust {} vs golden {}",
        windowed.len(),
        n_golden
    );

    let mut first_mismatch: Option<(usize, usize, i32, i32)> = None;
    let mut mismatch_count = 0usize;
    for i in 0..n_golden {
        let rust_w = &windowed[i].tensor;
        assert_eq!(rust_w.len(), 33 * 18, "tensor[{i}] length");
        for (k, &r) in rust_w.iter().enumerate() {
            let row = k / 18;
            let col = k % 18;
            let g = golden[(i, row, col)];
            if r != g {
                mismatch_count += 1;
                if first_mismatch.is_none() {
                    first_mismatch = Some((i, k, r, g));
                }
            }
        }
    }
    if let Some((i, k, r, g)) = first_mismatch {
        let row = k / 18;
        let col = k % 18;
        let pos_info = &windowed[i].pos_info;
        eprintln!(
            "[m7] first mismatch at cand {i} ({pos_info}) row={row} col={col}: rust={r} golden={g}"
        );
        eprintln!("[m7] total mismatches: {mismatch_count}");
    }
    assert_eq!(
        mismatch_count, 0,
        "M7 windowed tensor mismatches: {mismatch_count}"
    );
    eprintln!(
        "[m7] windowed tensor byte-identical: {} × 33 × 18 int32",
        n_golden
    );
}

#[test]
fn m7_pos_alt_tsv_byte_match() {
    let Some(dir) = golden_dir() else {
        eprintln!("[skip] CLAIR3_CONNECT_GOLDEN_DIR not set");
        return;
    };
    let tsv_path = dir.join("pileup_windowed_pos_alt.tsv");
    if !tsv_path.exists() {
        eprintln!("[skip] {} not present", tsv_path.display());
        return;
    }
    let bam = fixture_bam();
    let ref_ = fixture_ref();
    if !bam.exists() || !ref_.exists() {
        eprintln!("[skip] missing fixture BAM/ref");
        return;
    }

    let region = Region::parse(FIXTURE_REGION).expect("parse region");
    let params = PileupParams::default();
    let res = calculate_pileup(&region, &bam, &ref_, &params).expect("calculate_pileup");
    let windowed = window_candidates(&res.chunks, &res.candidates);

    let golden_text = std::fs::read_to_string(&tsv_path).unwrap();
    let golden_lines: Vec<&str> = golden_text.lines().collect();
    assert_eq!(
        windowed.len(),
        golden_lines.len(),
        "windowed line count: rust {} vs golden {}",
        windowed.len(),
        golden_lines.len()
    );

    fn lex_sort_i_entries(s: &str) -> String {

        let Some((depth, payload)) = s.split_once('-') else {
            return s.to_string();
        };
        let toks: Vec<&str> = payload.split(' ').filter(|t| !t.is_empty()).collect();
        let (mut xs, mut ds, mut is_, mut rs): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new());
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

    let mut mismatches = 0usize;
    for (i, wc) in windowed.iter().enumerate() {
        let rust_line = format!("{}\t{}", wc.pos_info, lex_sort_i_entries(&wc.alt_info));
        let parts: Vec<&str> = golden_lines[i].splitn(2, '\t').collect();
        let golden_normalized = if parts.len() == 2 {
            format!("{}\t{}", parts[0], lex_sort_i_entries(parts[1]))
        } else {
            golden_lines[i].to_string()
        };
        if rust_line != golden_normalized {
            mismatches += 1;
            if mismatches <= 5 {
                eprintln!("[m7-tsv] mismatch #{i}");
                eprintln!("         rust:   {rust_line}");
                eprintln!("         golden: {golden_normalized}");
            }
        }
    }
    assert_eq!(mismatches, 0, "M7 pos_alt TSV mismatches: {mismatches}");
    eprintln!(
        "[m7-tsv] all {} candidate lines byte-identical",
        windowed.len()
    );
}
