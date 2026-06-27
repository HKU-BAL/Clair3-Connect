use std::io::{BufRead, Write};

#[derive(Debug, Clone, Copy)]
pub struct MergeConfig {

    pub quality_score_for_pass: Option<f64>,

    pub print_ref: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MergeStats {

    pub fa_rows: usize,

    pub pileup_rows: usize,

    pub pileup_overridden: usize,
}

pub fn merge_vcf<RP: BufRead, RF: BufRead, W: Write>(
    pileup_vcf: RP,
    fa_vcf: RF,
    out: &mut W,
    ctg_filter: Option<&str>,
    config: MergeConfig,
) -> std::io::Result<MergeStats> {

    let mut fa_rows: Vec<(String, u64, String)> = Vec::new();
    let mut header_lines: Vec<String> = Vec::new();
    let mut buf = String::new();
    let mut fa_reader = fa_vcf;
    loop {
        buf.clear();
        if fa_reader.read_line(&mut buf)? == 0 {
            break;
        }
        let trimmed = buf.trim_end_matches(['\n', '\r']);
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('#') {
            header_lines.push(trimmed.to_string());
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 8 {
            continue;
        }
        if let Some(want) = ctg_filter {
            if cols[0] != want {
                continue;
            }
        }
        let pos: u64 = match cols[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let qual: f64 = match cols[5].parse() {
            Ok(q) => q,
            Err(_) => continue,
        };
        let ref_base = cols[3];
        let alt_base = cols[4];
        let is_reference = alt_base == "." || ref_base == alt_base;

        if is_reference && !config.print_ref {

            fa_rows.push((cols[0].to_string(), pos, String::new()));
            continue;
        }

        let row_text = if !is_reference {
            mark_low_qual(trimmed, config.quality_score_for_pass, qual)
        } else {
            trimmed.to_string()
        };
        fa_rows.push((cols[0].to_string(), pos, row_text));
    }

    let fa_pos_set: std::collections::HashSet<(String, u64)> = fa_rows
        .iter()
        .map(|(c, p, _)| (c.clone(), *p))
        .collect();

    let mut pileup_rows: Vec<(String, u64, String)> = Vec::new();
    let mut pileup_dropped: usize = 0;
    let mut pileup_reader = pileup_vcf;
    loop {
        buf.clear();
        if pileup_reader.read_line(&mut buf)? == 0 {
            break;
        }
        let trimmed = buf.trim_end_matches(['\n', '\r']);
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 8 {
            continue;
        }
        if let Some(want) = ctg_filter {
            if cols[0] != want {
                continue;
            }
        }
        let pos: u64 = match cols[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        if fa_pos_set.contains(&(cols[0].to_string(), pos)) {
            pileup_dropped += 1;
            continue;
        }
        let qual: f64 = match cols[5].parse() {
            Ok(q) => q,
            Err(_) => continue,
        };
        let ref_base = cols[3];
        let alt_base = cols[4];
        let is_reference = alt_base == "." || ref_base == alt_base;
        if is_reference && !config.print_ref {
            continue;
        }
        let row_text = if !is_reference {
            mark_low_qual(trimmed, config.quality_score_for_pass, qual)
        } else {
            trimmed.to_string()
        };
        pileup_rows.push((cols[0].to_string(), pos, row_text));
    }

    let mut merged: Vec<(String, u64, String)> = Vec::new();
    let mut fa_non_empty = 0;
    for (c, p, r) in fa_rows {
        if r.is_empty() {
            continue;
        }
        fa_non_empty += 1;
        merged.push((c, p, r));
    }
    let pileup_count = pileup_rows.len();
    merged.extend(pileup_rows);
    merged.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    for h in &header_lines {
        writeln!(out, "{}", h)?;
    }
    for (_, _, r) in &merged {
        writeln!(out, "{}", r)?;
    }

    Ok(MergeStats {
        fa_rows: fa_non_empty,
        pileup_rows: pileup_count,
        pileup_overridden: pileup_dropped,
    })
}

fn mark_low_qual(row: &str, quality_score_for_pass: Option<f64>, qual: f64) -> String {
    let Some(threshold) = quality_score_for_pass else {
        return row.to_string();
    };
    if qual > threshold {
        return row.to_string();
    }
    let mut cols: Vec<&str> = row.split('\t').collect();
    if cols.len() > 6 {
        cols[6] = "LowQual";
        cols.join("\t")
    } else {
        row.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pileup() -> &'static str {
        "\
##header
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
chr20\t100\t.\tA\tC\t10.0\tPASS\tP\tGT\t0/1
chr20\t200\t.\tA\tG\t20.0\tPASS\tP\tGT\t0/1
chr20\t300\t.\tA\tT\t30.0\tPASS\tP\tGT\t0/1
"
    }

    fn fa() -> &'static str {
        "\
##header
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
chr20\t200\t.\tA\tG\t50.0\tPASS\tF\tGT\t0/1
chr20\t400\t.\tA\tT\t60.0\tPASS\tF\tGT\t1/1
"
    }

    #[test]
    fn fa_overrides_pileup_at_same_position() {
        let mut out: Vec<u8> = Vec::new();
        let s = merge_vcf(
            pileup().as_bytes(),
            fa().as_bytes(),
            &mut out,
            Some("chr20"),
            MergeConfig {
                quality_score_for_pass: Some(2.0),
                print_ref: false,
            },
        )
        .unwrap();

        assert_eq!(s.fa_rows, 2);
        assert_eq!(s.pileup_rows, 2);
        assert_eq!(s.pileup_overridden, 1);
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();

        assert_eq!(lines.len(), 4);

        assert!(lines[1].contains("\tPASS\tF\t"));
        assert!(lines[1].contains("\t50.0\t"));
    }

    #[test]
    fn mark_low_qual_rewrites_filter() {
        let row = "chr20\t100\t.\tA\tC\t1.5\tPASS\tP\tGT\t0/1";
        let out = mark_low_qual(row, Some(2.0), 1.5);
        let cols: Vec<&str> = out.split('\t').collect();
        assert_eq!(cols[6], "LowQual");

        let row2 = "chr20\t100\t.\tA\tC\t10.0\tPASS\tP\tGT\t0/1";
        let out2 = mark_low_qual(row2, Some(2.0), 10.0);
        let cols2: Vec<&str> = out2.split('\t').collect();
        assert_eq!(cols2[6], "PASS");
    }

    #[test]
    fn print_ref_false_drops_ref_calls() {
        let pile = "\
#h
chr20\t100\t.\tA\t.\t50.0\tRefCall\t.\tGT\t0/0
chr20\t200\t.\tA\tC\t30.0\tPASS\t.\tGT\t0/1
";
        let fa = "#h\n";
        let mut out: Vec<u8> = Vec::new();
        let s = merge_vcf(
            pile.as_bytes(),
            fa.as_bytes(),
            &mut out,
            Some("chr20"),
            MergeConfig {
                quality_score_for_pass: None,
                print_ref: false,
            },
        )
        .unwrap();
        assert_eq!(s.pileup_rows, 1);
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(lines.len(), 1);
        assert!(lines[0].contains("\t200\t"));
    }
}
