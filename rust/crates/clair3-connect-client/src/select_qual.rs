use std::io::BufRead;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhaseQualSelection {

    pub qual_cut_off: f64,

    pub het_snp_count: usize,

    pub low_count: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QualSelection {

    pub var_qual_cut_off: f64,

    pub ref_qual_cut_off: f64,

    pub var_count: usize,

    pub ref_count: usize,
}

pub fn select_phase_qual<R: BufRead>(
    reader: R,
    low_qual_hete_var_pct: f64,
) -> std::io::Result<PhaseQualSelection> {
    let mut quals: Vec<f64> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() < 10 {
            continue;
        }
        let ref_base = cols[3];
        let alt_base = cols[4];
        if ref_base.len() != 1 || alt_base.len() != 1 {
            continue;
        }
        let qual = cols[5];
        let gt = cols[9].split(':').next().unwrap_or("");
        if gt != "0/1" {
            continue;
        }
        if let Ok(q) = qual.parse::<f64>() {
            quals.push(q);
        }
    }

    quals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = quals.len();

    let cutoff_idx = (low_qual_hete_var_pct * n as f64) as usize;
    let low_count = cutoff_idx.min(n);
    let qual_cut_off = if low_count == 0 {
        0.0
    } else {

        quals[low_count - 1]
    };
    Ok(PhaseQualSelection {
        qual_cut_off,
        het_snp_count: n,
        low_count,
    })
}

pub fn select_qual<R: BufRead>(
    reader: R,
    var_pct_full: f64,
    ref_pct_full: f64,
) -> std::io::Result<QualSelection> {
    let mut variant_qual_list: Vec<f64> = Vec::new();
    let mut ref_qual_list: Vec<f64> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() < 10 {
            continue;
        }
        let qual = cols[5];
        let gt = cols[9].split(':').next().unwrap_or("");
        let q = match qual.parse::<f64>() {
            Ok(q) => q,
            Err(_) => continue,
        };
        if gt == "0/0" {
            ref_qual_list.push(q);
        } else {
            variant_qual_list.push(q);
        }
    }

    ref_qual_list.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    variant_qual_list.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let var_cutoff_idx = (var_pct_full * variant_qual_list.len() as f64) as usize;
    let var_qual_cut_off = if var_cutoff_idx == 0 {
        0.0
    } else {
        variant_qual_list[var_cutoff_idx.min(variant_qual_list.len()) - 1]
    };

    let ref_cutoff_idx = (ref_pct_full * ref_qual_list.len() as f64) as usize;
    let ref_qual_cut_off = if ref_cutoff_idx == 0 {
        0.0
    } else {
        ref_qual_list[ref_cutoff_idx.min(ref_qual_list.len()) - 1]
    };

    Ok(QualSelection {
        var_qual_cut_off,
        ref_qual_cut_off,
        var_count: variant_qual_list.len(),
        ref_count: ref_qual_list.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_pileup_vcf() -> &'static str {
        "\
##fileformat=VCFv4.2
##source=Clair3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG003
chr20\t100\t.\tA\tC\t10.00\tPASS\tP\tGT:GQ\t0/1:10
chr20\t200\t.\tA\tG\t20.00\tPASS\tP\tGT:GQ\t0/1:20
chr20\t300\t.\tA\tT\t30.00\tPASS\tP\tGT:GQ\t0/1:30
chr20\t400\t.\tA\tC\t40.00\tPASS\tP\tGT:GQ\t1/1:40
chr20\t500\t.\tA\tC\t50.00\tPASS\tP\tGT:GQ\t0/1:50
chr20\t600\t.\tA\t.\t5.00\tRefCall\tP\tGT:GQ\t0/0:5
chr20\t700\t.\tA\tACG\t30.00\tPASS\tP\tGT:GQ\t0/1:30
"
    }

    #[test]
    fn select_phase_qual_excludes_indels_and_non_het() {
        let bytes = sample_pileup_vcf().as_bytes();
        let s = select_phase_qual(bytes, 0.3).unwrap();

        assert_eq!(s.het_snp_count, 4);
        assert_eq!(s.low_count, 1);
        assert!((s.qual_cut_off - 10.0).abs() < 1e-9);
    }

    #[test]
    fn select_phase_qual_zero_when_no_het_snps() {
        let bytes = b"#header\nchr20\t100\t.\tA\tC\t10\tPASS\t.\tGT\t0/0\n".as_slice();
        let s = select_phase_qual(bytes, 0.3).unwrap();
        assert_eq!(s.het_snp_count, 0);
        assert_eq!(s.qual_cut_off, 0.0);
    }

    #[test]
    fn select_qual_separates_ref_and_variant() {
        let bytes = sample_pileup_vcf().as_bytes();
        let s = select_qual(bytes, 0.3, 0.3).unwrap();

        assert_eq!(s.var_count, 6);
        assert!((s.var_qual_cut_off - 10.0).abs() < 1e-9);

        assert_eq!(s.ref_count, 1);
        assert_eq!(s.ref_qual_cut_off, 0.0);
    }

    #[test]
    fn select_phase_qual_index_floor() {

        let mut s = String::from("##h\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n");
        for q in 1..=10 {
            s.push_str(&format!("chr20\t{}\t.\tA\tC\t{}.0\tPASS\t.\tGT\t0/1\n", q, q));
        }

        let r = select_phase_qual(s.as_bytes(), 0.3).unwrap();
        assert_eq!(r.low_count, 3);
        assert!((r.qual_cut_off - 3.0).abs() < 1e-9);

        let r0 = select_phase_qual(s.as_bytes(), 0.0).unwrap();
        assert_eq!(r0.qual_cut_off, 0.0);

        let r1 = select_phase_qual(s.as_bytes(), 1.0).unwrap();
        assert!((r1.qual_cut_off - 10.0).abs() < 1e-9);
    }
}
