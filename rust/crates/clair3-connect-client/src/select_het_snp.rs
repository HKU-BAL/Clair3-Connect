use std::io::{BufRead, Write};

#[derive(Debug, Clone, Copy)]
pub enum HetSnpFilterMode {

    QualCutoff(f64),

    BottomPercent(f64),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FilterStats {

    pub kept: usize,

    pub dropped: usize,

    pub total_het_snp: usize,
}

pub fn filter_het_snp_phasing<R: BufRead, W: Write>(
    mut input: R,
    output: &mut W,
    ctg_filter: Option<&str>,
    mode: HetSnpFilterMode,
) -> std::io::Result<FilterStats> {

    let mut header_lines: Vec<String> = Vec::new();

    let mut entries: Vec<(u64, String, f64)> = Vec::new();

    let mut buf = String::new();
    loop {
        buf.clear();
        let n = input.read_line(&mut buf)?;
        if n == 0 {
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
        if cols.len() < 10 {
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
        let ref_base = cols[3];
        let alt_base = cols[4];
        if ref_base.len() != 1 || alt_base.len() != 1 {
            continue;
        }

        let gt_raw = cols[9].split(':').next().unwrap_or("");
        let gt_norm: String = gt_raw.chars().map(|c| if c == '|' { '/' } else { c }).collect();
        if gt_norm != "0/1" && gt_norm != "1/0" {
            continue;
        }
        let qual: f64 = match cols[5].parse() {
            Ok(q) => q,
            Err(_) => continue,
        };

        let mut owned_cols: Vec<String> = cols.iter().map(|s| s.to_string()).collect();
        let fmt_idx = owned_cols.len() - 2;
        let sample_idx = owned_cols.len() - 1;
        if owned_cols[fmt_idx].contains("PL") {
            let fmt_parts: Vec<&str> = owned_cols[fmt_idx].split(':').collect();
            if let Some(pl_pos) = fmt_parts.iter().position(|f| *f == "PL") {
                let new_fmt: Vec<&str> = fmt_parts
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != pl_pos)
                    .map(|(_, s)| *s)
                    .collect();
                let new_fmt_str = new_fmt.join(":");

                let sample_parts: Vec<&str> = owned_cols[sample_idx].split(':').collect();
                let new_sample: Vec<&str> = sample_parts
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != pl_pos)
                    .map(|(_, s)| *s)
                    .collect();
                let new_sample_str = new_sample.join(":");

                owned_cols[fmt_idx] = new_fmt_str;
                owned_cols[sample_idx] = new_sample_str;
            }
        }
        let row = owned_cols.join("\t");
        entries.push((pos, row, qual));
    }

    let total_het_snp = entries.len();
    let drop_positions: std::collections::HashSet<u64> = match mode {
        HetSnpFilterMode::QualCutoff(cutoff) => entries
            .iter()
            .filter(|(_, _, q)| *q < cutoff)
            .map(|(p, _, _)| *p)
            .collect(),
        HetSnpFilterMode::BottomPercent(pct) => {
            let mut sorted: Vec<(u64, f64)> = entries.iter().map(|(p, _, q)| (*p, *q)).collect();
            sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            let drop_n = (pct * sorted.len() as f64) as usize;
            sorted.into_iter().take(drop_n).map(|(p, _)| p).collect()
        }
    };

    let kept_rows: Vec<&(u64, String, f64)> = {
        let mut v: Vec<_> = entries.iter().filter(|(p, _, _)| !drop_positions.contains(p)).collect();
        v.sort_by_key(|(p, _, _)| *p);
        v
    };
    let kept = kept_rows.len();
    let dropped = total_het_snp - kept;

    for h in &header_lines {
        writeln!(output, "{}", h)?;
    }
    for (_, row, _) in &kept_rows {
        writeln!(output, "{}", row)?;
    }

    Ok(FilterStats {
        kept,
        dropped,
        total_het_snp,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vcf_text() -> &'static str {
        "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG003
chr20\t100\t.\tA\tC\t10.00\tPASS\tP\tGT:GQ:PL\t0/1:10:100,0,200
chr20\t300\t.\tA\tT\t30.00\tPASS\tP\tGT:GQ:PL\t0/1:30:100,0,200
chr20\t200\t.\tA\tG\t20.00\tPASS\tP\tGT:GQ\t1|0:20
chr20\t400\t.\tA\tC\t40.00\tPASS\tP\tGT:GQ:PL\t1/1:40:200,0,0
chr20\t500\t.\tA\tC\t50.00\tPASS\tP\tGT:GQ:PL\t0/1:50:300,0,100
"
    }

    #[test]
    fn qual_cutoff_drops_below_threshold() {
        let mut out: Vec<u8> = Vec::new();
        let s = filter_het_snp_phasing(
            vcf_text().as_bytes(),
            &mut out,
            Some("chr20"),
            HetSnpFilterMode::QualCutoff(20.0),
        )
        .unwrap();

        assert_eq!(s.total_het_snp, 4);
        assert_eq!(s.kept, 3);
        let text = String::from_utf8(out).unwrap();

        let lines: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(lines.len(), 3);
        assert!(lines[0].starts_with("chr20\t200"));
        assert!(lines[1].starts_with("chr20\t300"));
        assert!(lines[2].starts_with("chr20\t500"));
    }

    #[test]
    fn pl_field_stripped() {
        let mut out: Vec<u8> = Vec::new();
        filter_het_snp_phasing(
            vcf_text().as_bytes(),
            &mut out,
            Some("chr20"),
            HetSnpFilterMode::QualCutoff(0.0),
        )
        .unwrap();
        let text = String::from_utf8(out).unwrap();

        for line in text.lines() {
            if line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            assert!(!cols[8].contains("PL"), "FORMAT still has PL: {line}");

        }
    }

    #[test]
    fn collapses_phased_to_unphased() {

        let mut out: Vec<u8> = Vec::new();
        let s = filter_het_snp_phasing(
            vcf_text().as_bytes(),
            &mut out,
            None,
            HetSnpFilterMode::QualCutoff(0.0),
        )
        .unwrap();

        assert_eq!(s.total_het_snp, 4);
    }

    #[test]
    fn bottom_percent_mode_drops_lowest() {
        let mut out: Vec<u8> = Vec::new();
        let s = filter_het_snp_phasing(
            vcf_text().as_bytes(),
            &mut out,
            None,
            HetSnpFilterMode::BottomPercent(0.5),
        )
        .unwrap();

        assert_eq!(s.kept, 2);
    }
}
