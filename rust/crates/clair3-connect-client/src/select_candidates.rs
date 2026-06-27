use std::io::BufRead;

#[derive(Debug, Clone, Copy)]
pub enum CandidateFilterMode {

    QualCutoffs {

        var_cutoff: f64,

        ref_cutoff: f64,
    },

    BottomPercent {

        var_pct_full: f64,

        ref_pct_full: f64,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct CandidateSelection {

    pub positions: Vec<u64>,

    pub variants_total: usize,

    pub variants_kept: usize,

    pub refs_total: usize,

    pub refs_kept: usize,
}

pub fn select_candidates<R: BufRead>(
    reader: R,
    ctg_filter: Option<&str>,
    mode: CandidateFilterMode,
) -> std::io::Result<CandidateSelection> {
    let mut ref_call_pos_list: Vec<(u64, f64)> = Vec::new();
    let mut variant_pos_list: Vec<(u64, f64)> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 6 {
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
        let qual: f64 = match cols[5].parse() {
            Ok(q) => q,
            Err(_) => continue,
        };

        if alt_base == "." || alt_base == ref_base {
            ref_call_pos_list.push((pos, qual));
        } else {
            variant_pos_list.push((pos, qual));
        }
    }

    let refs_total = ref_call_pos_list.len();
    let variants_total = variant_pos_list.len();

    let (kept_refs, kept_vars) = match mode {
        CandidateFilterMode::QualCutoffs { var_cutoff, ref_cutoff } => {
            let kept_refs: Vec<u64> = ref_call_pos_list
                .iter()
                .filter(|(_, q)| *q < ref_cutoff)
                .map(|(p, _)| *p)
                .collect();
            let kept_vars: Vec<u64> = variant_pos_list
                .iter()
                .filter(|(_, q)| *q < var_cutoff)
                .map(|(p, _)| *p)
                .collect();
            (kept_refs, kept_vars)
        }
        CandidateFilterMode::BottomPercent { var_pct_full, ref_pct_full } => {
            let mut refs_sorted = ref_call_pos_list.clone();
            refs_sorted.sort_by(|a, b| {
                a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
            });
            let ref_n = (ref_pct_full * refs_sorted.len() as f64) as usize;
            let kept_refs: Vec<u64> =
                refs_sorted.into_iter().take(ref_n).map(|(p, _)| p).collect();

            let mut vars_sorted = variant_pos_list.clone();
            vars_sorted.sort_by(|a, b| {
                a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
            });
            let var_n = (var_pct_full * vars_sorted.len() as f64) as usize;
            let kept_vars: Vec<u64> =
                vars_sorted.into_iter().take(var_n).map(|(p, _)| p).collect();

            (kept_refs, kept_vars)
        }
    };

    let refs_kept = kept_refs.len();
    let variants_kept = kept_vars.len();

    let mut all: std::collections::BTreeSet<u64> = std::collections::BTreeSet::new();
    for p in kept_refs {
        all.insert(p);
    }
    for p in kept_vars {
        all.insert(p);
    }

    Ok(CandidateSelection {
        positions: all.into_iter().collect(),
        variants_total,
        variants_kept,
        refs_total,
        refs_kept,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vcf_text() -> &'static str {
        "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG003
chr20\t100\t.\tA\tC\t10.00\tPASS\tP\tGT:GQ\t0/1:10
chr20\t200\t.\tA\tG\t20.00\tPASS\tP\tGT:GQ\t0/1:20
chr20\t300\t.\tA\t.\t5.00\tRefCall\tP\tGT:GQ\t0/0:5
chr20\t400\t.\tA\t.\t30.00\tRefCall\tP\tGT:GQ\t0/0:30
chr20\t500\t.\tA\tT\t50.00\tPASS\tP\tGT:GQ\t1/1:50
chr20\t600\t.\tA\tA\t0.00\tRefCall\tP\tGT:GQ\t0/0:0
chr1\t900\t.\tA\tC\t99.00\tPASS\tP\tGT:GQ\t0/1:99
"
    }

    #[test]
    fn quality_cutoffs_filter_below_threshold() {
        let s = select_candidates(
            vcf_text().as_bytes(),
            Some("chr20"),
            CandidateFilterMode::QualCutoffs {
                var_cutoff: 30.0,
                ref_cutoff: 10.0,
            },
        )
        .unwrap();

        assert_eq!(s.positions, vec![100, 200, 300, 600]);
        assert_eq!(s.variants_kept, 2);
        assert_eq!(s.refs_kept, 2);

        assert_eq!(s.variants_total, 3);
        assert_eq!(s.refs_total, 3);
    }

    #[test]
    fn bottom_percent_mode_picks_lowest_qual() {
        let s = select_candidates(
            vcf_text().as_bytes(),
            Some("chr20"),
            CandidateFilterMode::BottomPercent {
                var_pct_full: 0.34,
                ref_pct_full: 0.66,
            },
        )
        .unwrap();

        assert_eq!(s.positions, vec![100, 600]);
    }

    #[test]
    fn refcall_dot_alt_treated_as_ref() {
        let bytes = "#h\nchr20\t100\t.\tA\t.\t5\tRefCall\t.\tGT\t0/0\n".as_bytes();
        let s = select_candidates(
            bytes,
            None,
            CandidateFilterMode::QualCutoffs { var_cutoff: 100.0, ref_cutoff: 100.0 },
        )
        .unwrap();
        assert_eq!(s.refs_total, 1);
        assert_eq!(s.variants_total, 0);
        assert_eq!(s.positions, vec![100]);
    }

    #[test]
    fn dedup_when_pos_in_both_lists() {

        let bytes = "\
#h
chr20\t100\t.\tA\tC\t5\tPASS\t.\tGT\t0/1
chr20\t100\t.\tA\t.\t3\tRefCall\t.\tGT\t0/0
"
        .as_bytes();
        let s = select_candidates(
            bytes,
            None,
            CandidateFilterMode::QualCutoffs { var_cutoff: 100.0, ref_cutoff: 100.0 },
        )
        .unwrap();
        assert_eq!(s.positions, vec![100]);
    }
}
