use clair3_connect_core::full_alignment::driver::{calculate_full_alignment_with_haplotag, FaInputs};
use clair3_connect_core::full_alignment::phase_block_index::Variant;
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
    assert!(header.contains("'|i1'") || header.contains("'int8'"));
    bytes[header_len..].iter().map(|b| *b as i8).collect()
}

fn read_variants_tsv(path: &std::path::Path) -> Vec<Variant> {
    std::fs::read_to_string(path)
        .expect("read variants tsv")
        .lines()
        .map(|line| {
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(cols.len(), 5, "expected 5 columns, got {cols:?}");
            Variant {
                position: cols[0].parse().expect("position"),
                ref_base: cols[1].as_bytes()[0],
                alt_base: cols[2].as_bytes()[0],
                genotype: cols[3].parse().expect("genotype"),
                phase_set: cols[4].parse().expect("phase_set"),
            }
        })
        .collect()
}

fn read_candidates(path: &std::path::Path) -> Vec<u64> {
    std::fs::read_to_string(path)
        .expect("read candidates")
        .lines()
        .map(|l| l.parse().expect("u64"))
        .collect()
}

#[test]
fn fa_haplotag_matrix_bytes_match_libclair3() {
    let arts = artifacts_dir();
    let cand_path = arts.join("fa_haplotag_candidates.tsv");
    let matrix_path = arts.join("fa_haplotag_matrix.npy");
    let variants_path = arts.join("fa_haplotag_variants.tsv");
    if !cand_path.exists() || !matrix_path.exists() || !variants_path.exists() {
        eprintln!("Skipping — Phase 4b fixture not present");
        return;
    }
    let candidates = read_candidates(&cand_path);
    let variants = read_variants_tsv(&variants_path);
    eprintln!("loaded {} candidates, {} variants", candidates.len(), variants.len());

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
    let result = calculate_full_alignment_with_haplotag(&inputs, &variants)
        .expect("calculate_full_alignment_with_haplotag");

    let expected = read_npy_i8(&matrix_path);
    assert_eq!(result.matrix.len(), expected.len(), "matrix size differs");

    let mut mismatches = 0;
    let mut by_channel = [0_usize; 8];
    let mut first_diff = None;
    for (i, (a, b)) in result.matrix.iter().zip(expected.iter()).enumerate() {
        if a != b {
            mismatches += 1;
            let ch = i % 8;
            by_channel[ch] += 1;
            if first_diff.is_none() {
                first_diff = Some((i, *a, *b));
            }
        }
    }
    eprintln!("haplotag matrix mismatches total: {mismatches}");
    eprintln!("by_channel: {:?} (idx: ref alt strand mq bq af ins hap)", by_channel);
    if let Some((i, got, exp)) = first_diff {
        let off_can = 89 * 33 * 8;
        let off_row = 33 * 8;
        let off_col = 8;
        let c = i / off_can;
        let d = (i % off_can) / off_row;
        let p = (i % off_row) / off_col;
        let ch = i % off_col;
        eprintln!(
            "first diff at idx {i} = candidate={c} depth_row={d} pos_in_window={p} channel={ch}: \
             got={got} exp={exp}"
        );
    }
    assert_eq!(mismatches, 0, "Phase 4b G6 byte-equal failed");
}

#[test]
fn fa_haplotag_alt_info_match_libclair3() {
    let arts = artifacts_dir();
    let alt_info_path = arts.join("fa_haplotag_alt_info.tsv");
    let cand_path = arts.join("fa_haplotag_candidates.tsv");
    let variants_path = arts.join("fa_haplotag_variants.tsv");
    if !alt_info_path.exists() || !cand_path.exists() || !variants_path.exists() {
        eprintln!("Skipping — Phase 4b fixture not present");
        return;
    }
    let candidates = read_candidates(&cand_path);
    let variants = read_variants_tsv(&variants_path);
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
    let result = calculate_full_alignment_with_haplotag(&inputs, &variants)
        .expect("calculate_full_alignment_with_haplotag");

    let expected: Vec<String> = std::fs::read_to_string(&alt_info_path)
        .unwrap()
        .lines()
        .map(|l| l.to_string())
        .collect();
    assert_eq!(result.alt_info_strings.len(), expected.len());

    let mut mismatches = 0;
    let mut first = None;
    for (i, (got, exp)) in result.alt_info_strings.iter().zip(expected.iter()).enumerate() {
        let g = got.trim_end();
        let e = exp.trim_end();
        if g != e {
            mismatches += 1;
            if first.is_none() {
                first = Some((i, g.to_string(), e.to_string()));
            }
        }
    }
    if let Some((i, g, e)) = &first {
        eprintln!("first alt_info diff row {i}:\n  got: {g:?}\n  exp: {e:?}");
    }
    assert_eq!(mismatches, 0, "Phase 4b G8 alt_info parity failed");
}
