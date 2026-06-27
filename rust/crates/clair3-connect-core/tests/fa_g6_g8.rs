use clair3_connect_core::full_alignment::driver::{calculate_full_alignment_no_haplotag, FaInputs};
use std::path::PathBuf;

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("test_data/ont_chr20_demo")
}

fn artifacts_dir() -> PathBuf {
    fixture_dir().join("artifacts")
}

fn read_npy_i8(path: &std::path::Path) -> Vec<i8> {
    let bytes = std::fs::read(path).expect("read npy");
    assert_eq!(&bytes[0..6], b"\x93NUMPY");
    let major = bytes[6];
    let header_len = if major == 1 {
        u16::from_le_bytes([bytes[8], bytes[9]]) as usize + 10
    } else {
        u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize + 12
    };
    let header = std::str::from_utf8(&bytes[10..header_len]).unwrap();
    assert!(
        header.contains("'|i1'") || header.contains("'int8'"),
        "expected int8 npy: {header}"
    );
    bytes[header_len..].iter().map(|b| *b as i8).collect()
}

#[test]
fn fa_no_haplotag_matrix_bytes_match_libclair3() {
    let arts = artifacts_dir();
    let cand_path = arts.join("fa_no_haplotag_candidates.tsv");
    let matrix_path = arts.join("fa_no_haplotag_matrix.npy");
    if !cand_path.exists() || !matrix_path.exists() {
        eprintln!("Skipping — Phase 0 FA fixture not present");
        return;
    }
    let candidates: Vec<u64> = std::fs::read_to_string(&cand_path)
        .expect("read candidates")
        .lines()
        .map(|l| l.parse().unwrap())
        .collect();
    assert!(!candidates.is_empty());

    let bam_path = fixture_dir().join("HG003_chr20_demo.bam");
    let fa_path = fixture_dir().join("GRCh38_no_alt_chr20.fa");

    let inputs = FaInputs {
        bam_path: &bam_path,
        fasta_path: &fa_path,

        region: ("chr20".to_string(), 99999, 300000),
        candidates: &candidates,
        min_mq: 5,
        matrix_depth: 89,
        max_indel_length: 50,
    };
    let result =
        calculate_full_alignment_no_haplotag(&inputs).expect("calculate_full_alignment_no_haplotag");

    let expected_matrix = read_npy_i8(&matrix_path);
    eprintln!(
        "rust_matrix.len={} expected_matrix.len={} candidate_num={}",
        result.matrix.len(),
        expected_matrix.len(),
        result.candidate_num
    );
    assert_eq!(result.matrix.len(), expected_matrix.len(), "matrix size differs");

    let mut mismatches = 0;
    let mut first_diff = None;
    for (i, (a, b)) in result.matrix.iter().zip(expected_matrix.iter()).enumerate() {
        if a != b {
            mismatches += 1;
            if first_diff.is_none() {
                first_diff = Some((i, *a, *b));
            }
        }
    }
    eprintln!("matrix mismatches: {mismatches}");
    if let Some((i, got, exp)) = first_diff {

        let off_can = result.matrix_depth * result.no_of_positions * result.channel_size;
        let off_row = result.no_of_positions * result.channel_size;
        let off_col = result.channel_size;
        let c = i / off_can;
        let d = (i % off_can) / off_row;
        let p = (i % off_row) / off_col;
        let ch = i % off_col;
        eprintln!(
            "first diff at idx {i} = candidate={c} depth_row={d} pos_in_window={p} channel={ch}: \
             got={got} exp={exp}",
        );
    }
    assert_eq!(mismatches, 0, "G6 byte-equal failed");
}

#[test]
fn fa_no_haplotag_alt_info_match_libclair3() {
    let arts = artifacts_dir();
    let cand_path = arts.join("fa_no_haplotag_candidates.tsv");
    let alt_info_path = arts.join("fa_no_haplotag_alt_info.tsv");
    if !cand_path.exists() || !alt_info_path.exists() {
        eprintln!("Skipping — Phase 0 FA fixture not present");
        return;
    }
    let candidates: Vec<u64> = std::fs::read_to_string(&cand_path)
        .expect("read candidates")
        .lines()
        .map(|l| l.parse().unwrap())
        .collect();
    let bam_path = fixture_dir().join("HG003_chr20_demo.bam");
    let fa_path = fixture_dir().join("GRCh38_no_alt_chr20.fa");

    let inputs = FaInputs {
        bam_path: &bam_path,
        fasta_path: &fa_path,
        region: ("chr20".to_string(), 99999, 300000),
        candidates: &candidates,
        min_mq: 5,
        matrix_depth: 89,
        max_indel_length: 50,
    };
    let result =
        calculate_full_alignment_no_haplotag(&inputs).expect("calculate_full_alignment_no_haplotag");

    let expected_lines: Vec<String> = std::fs::read_to_string(&alt_info_path)
        .expect("read alt_info")
        .lines()
        .map(|l| l.to_string())
        .collect();
    assert_eq!(result.alt_info_strings.len(), expected_lines.len(),
               "alt_info row count differs");

    let mut mismatches = 0;
    let mut first: Option<(usize, String, String)> = None;
    for (i, (got, exp)) in result
        .alt_info_strings
        .iter()
        .zip(expected_lines.iter())
        .enumerate()
    {

        let got_t = got.trim_end();
        let exp_t = exp.trim_end();
        if got_t != exp_t {
            mismatches += 1;
            if first.is_none() {
                first = Some((i, got_t.to_string(), exp_t.to_string()));
            }
        }
    }
    if let Some((i, g, e)) = &first {
        eprintln!("first alt_info diff at row {i}:\n  got: {g:?}\n  exp: {e:?}");
    }
    assert_eq!(mismatches, 0, "G8 alt_info parity failed");
}
