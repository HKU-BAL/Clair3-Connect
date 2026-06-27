use clair3_connect_core::vcf::find_alt_base::decode_alt_info;
use clair3_connect_core::vcf::output_from::output_from;
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
fn pileup_output_from_matches_python() {
    let dir = fixture_dir();
    let y_path = dir.join("Y_pileup_cpu.npy");
    let pos_alt_path = dir.join("pileup_windowed_pos_alt.tsv");
    let parity_path = dir.join("pileup_output_from.tsv");

    if !y_path.exists() || !parity_path.exists() {
        eprintln!("Skipping — Phase 0 fixture not present at {:?}", dir);
        return;
    }

    let (n_rows, y_flat) = read_npy_f32(&y_path, Y_COLS);
    let pos_alt_text = std::fs::read_to_string(&pos_alt_path).expect("read pos_alt");
    let pos_alt_lines: Vec<&str> = pos_alt_text.lines().collect();
    let parity_text = std::fs::read_to_string(&parity_path).expect("read parity");
    let parity_lines: Vec<&str> = parity_text.lines().collect();
    assert_eq!(n_rows, pos_alt_lines.len());
    assert_eq!(n_rows, parity_lines.len());

    let mut total = 0;
    let mut by_class: std::collections::BTreeMap<String, usize> =
        std::collections::BTreeMap::new();
    let mut mismatches = 0;
    let mut first_mismatch_msg: Option<String> = None;

    for row in 0..n_rows {
        let gt21 = &y_flat[row * Y_COLS..row * Y_COLS + 21];
        let geno = &y_flat[row * Y_COLS + 21..row * Y_COLS + 24];
        let pos_alt_line = pos_alt_lines[row];
        let (chr_pos_ref, payload) = pos_alt_line.split_once('\t').unwrap();
        let ref_field = chr_pos_ref.rsplit(':').next().unwrap();
        let (_depth, entries) = decode_alt_info(payload).expect("decode_alt_info");

        let d = output_from(ref_field, 0, gt21, geno, &entries, 50);

        let rust_flags = format!(
            "{}{}{}{}{}{}{}{}{}{}",
            u8::from(d.is_reference),
            u8::from(d.is_homo_snp),
            u8::from(d.is_hetero_snp),
            u8::from(d.is_homo_insertion),
            u8::from(d.is_hetero_acgt_ins),
            u8::from(d.is_hetero_insins),
            u8::from(d.is_homo_deletion),
            u8::from(d.is_hetero_acgt_del),
            u8::from(d.is_hetero_deldel),
            u8::from(d.is_insertion_and_deletion),
        );

        let parity_line = parity_lines[row];
        let parts: Vec<&str> = parity_line.split('\t').collect();
        assert_eq!(parts.len(), 4, "row {row}: bad parity line {parity_line:?}");
        let exp_flags = parts[0];
        let exp_ref = parts[1];
        let exp_alt = parts[2];
        let exp_p_hex = parts[3];

        let bytes: Vec<u8> = (0..exp_p_hex.len() / 2)
            .map(|i| u8::from_str_radix(&exp_p_hex[i * 2..i * 2 + 2], 16).unwrap())
            .collect();
        let exp_p =
            f32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);

        let class = {
            let mut s = String::new();
            for (i, c) in exp_flags.chars().enumerate() {
                if c == '1' {
                    s.push_str(match i {
                        0 => "ref",
                        1 => "homoSNP",
                        2 => "hetSNP",
                        3 => "homoIns",
                        4 => "hetACGTIns",
                        5 => "hetInsIns",
                        6 => "homoDel",
                        7 => "hetACGTDel",
                        8 => "hetDelDel",
                        9 => "insDel",
                        _ => unreachable!(),
                    });
                    break;
                }
            }
            if s.is_empty() {
                "noClass".to_string()
            } else {
                s
            }
        };
        *by_class.entry(class.clone()).or_insert(0) += 1;
        total += 1;

        let p_bits_match = d.maximum_probability.to_bits() == exp_p.to_bits();
        let ok = rust_flags == exp_flags
            && d.reference_base == exp_ref
            && d.alternate_base == exp_alt
            && p_bits_match;
        if !ok {
            mismatches += 1;
            if first_mismatch_msg.is_none() {
                first_mismatch_msg = Some(format!(
                    "row {row} class={class}\n  exp: flags={exp_flags} ref={exp_ref} alt={exp_alt} p={exp_p:e} bits={:#x}\n  got: flags={rust_flags} ref={} alt={} p={:e} bits={:#x}",
                    exp_p.to_bits(),
                    d.reference_base,
                    d.alternate_base,
                    d.maximum_probability,
                    d.maximum_probability.to_bits(),
                ));
            }
        }
    }

    eprintln!("total={total} mismatches={mismatches}");
    eprintln!("by_class: {by_class:?}");
    if mismatches > 0 {
        panic!("output_from parity failures: {mismatches}\nfirst:\n{}",
               first_mismatch_msg.unwrap());
    }
}
