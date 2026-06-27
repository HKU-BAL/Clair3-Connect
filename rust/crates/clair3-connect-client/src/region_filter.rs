use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Interval0Based {

    pub start: u64,

    pub end: u64,
}

impl Interval0Based {

    pub fn new(start: u64, end: u64) -> Self {
        debug_assert!(start < end, "empty interval [{start}, {end})");
        Self { start, end }
    }

    pub fn len(&self) -> u64 {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.end <= self.start
    }

    pub fn contains(&self, pos: u64) -> bool {
        pos >= self.start && pos < self.end
    }

    pub fn intersect_range(&self, s: u64, e: u64) -> Option<Interval0Based> {
        let start = self.start.max(s);
        let end = self.end.min(e);
        if start < end {
            Some(Interval0Based { start, end })
        } else {
            None
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct RegionFilter {

    by_contig: HashMap<String, Vec<Interval0Based>>,
}

impl RegionFilter {

    pub fn empty() -> Self {
        Self {
            by_contig: HashMap::new(),
        }
    }

    pub fn is_unrestricted(&self) -> bool {
        self.by_contig.is_empty()
    }

    pub fn contigs(&self) -> impl Iterator<Item = &str> {
        self.by_contig.keys().map(|s| s.as_str())
    }

    pub fn covered_bp(&self, contig: &str) -> u64 {
        self.by_contig
            .get(contig)
            .map(|ivs| ivs.iter().map(|i| i.len()).sum())
            .unwrap_or(0)
    }

    pub fn intervals_for(&self, contig: &str) -> Option<&[Interval0Based]> {
        self.by_contig.get(contig).map(|v| v.as_slice())
    }

    fn push_raw(&mut self, contig: &str, start: u64, end: u64) {
        if start >= end {
            return;
        }
        self.by_contig
            .entry(contig.to_string())
            .or_default()
            .push(Interval0Based::new(start, end));
    }

    fn finalize(&mut self) {
        for ivs in self.by_contig.values_mut() {
            if ivs.is_empty() {
                continue;
            }
            ivs.sort();
            let mut merged: Vec<Interval0Based> = Vec::with_capacity(ivs.len());
            for iv in ivs.iter() {
                if let Some(last) = merged.last_mut() {
                    if iv.start <= last.end {
                        last.end = last.end.max(iv.end);
                        continue;
                    }
                }
                merged.push(*iv);
            }
            *ivs = merged;
        }
    }

    pub fn from_bed_path(path: &Path) -> std::io::Result<Self> {
        let text = std::fs::read_to_string(path)?;
        Self::from_bed_text(&text)
    }

    pub fn from_bed_text(text: &str) -> std::io::Result<Self> {
        let mut f = Self::empty();
        for (lineno, line) in text.lines().enumerate() {
            let trimmed = line.trim_start();
            if trimmed.is_empty()
                || trimmed.starts_with('#')
                || trimmed.starts_with("track")
                || trimmed.starts_with("browser")
            {
                continue;
            }
            let cols: Vec<&str> = trimmed.split('\t').collect();
            if cols.len() < 3 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "bed line {} has {} columns; need >= 3 (chrom, start, end)",
                        lineno + 1,
                        cols.len()
                    ),
                ));
            }
            let chrom = cols[0];
            let start: u64 = cols[1].parse().map_err(|e| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("bed line {} start: {e}", lineno + 1),
                )
            })?;
            let end: u64 = cols[2].parse().map_err(|e| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("bed line {} end: {e}", lineno + 1),
                )
            })?;
            f.push_raw(chrom, start, end);
        }
        f.finalize();
        Ok(f)
    }

    pub fn from_region_strings(
        specs: &[&str],
        contig_lengths: &HashMap<String, u64>,
    ) -> std::io::Result<Self> {
        let mut f = Self::empty();
        for spec in specs {
            let s = spec.trim();
            if s.is_empty() {
                continue;
            }
            if let Some((chrom, range)) = s.split_once(':') {
                let (a, b) = range.split_once('-').ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("region {s:?}: expected chrom:start-end"),
                    )
                })?;
                let start_1b: u64 = a.replace(',', "").parse().map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("region {s:?} start: {e}"),
                    )
                })?;
                let end_1b: u64 = b.replace(',', "").parse().map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("region {s:?} end: {e}"),
                    )
                })?;
                if start_1b == 0 || end_1b == 0 {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "region {s:?}: 1-based coordinates must be >= 1"
                        ),
                    ));
                }
                if start_1b > end_1b {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("region {s:?}: start > end"),
                    ));
                }

                f.push_raw(chrom, start_1b - 1, end_1b);
            } else {
                let length = contig_lengths.get(s).copied().ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("region {s:?}: contig not in fai"),
                    )
                })?;
                f.push_raw(s, 0, length);
            }
        }
        f.finalize();
        Ok(f)
    }

    pub fn union(mut self, other: Self) -> Self {
        for (contig, ivs) in other.by_contig {
            let dst = self.by_contig.entry(contig).or_default();
            dst.extend(ivs);
        }
        self.finalize();
        self
    }

    pub fn contains_position(&self, contig: &str, pos_0b: u64) -> bool {
        let ivs = match self.by_contig.get(contig) {
            Some(v) => v,
            None => return false,
        };

        ivs.binary_search_by(|iv| {
            if iv.end <= pos_0b {
                std::cmp::Ordering::Less
            } else if iv.start > pos_0b {
                std::cmp::Ordering::Greater
            } else {
                std::cmp::Ordering::Equal
            }
        })
        .is_ok()
    }

    pub fn intersect_chunk_1based(
        &self,
        contig: &str,
        chunk_start_1b: u64,
        chunk_end_1b: u64,
    ) -> Vec<(u64, u64)> {

        let ivs = match self.by_contig.get(contig) {
            Some(v) => v.as_slice(),
            None => return Vec::new(),
        };
        let chunk_start_0b = chunk_start_1b.saturating_sub(1);
        let chunk_end_0b = chunk_end_1b;
        let mut out = Vec::new();
        for iv in ivs {

            if iv.start >= chunk_end_0b {
                break;
            }
            if iv.end <= chunk_start_0b {
                continue;
            }
            if let Some(int) = iv.intersect_range(chunk_start_0b, chunk_end_0b) {

                let s = int.start + 1;
                let e = int.end;
                out.push((s, e));
            }
        }
        out
    }
}

pub fn filter_vcf_by_region(buf: &[u8], filter: &RegionFilter) -> Vec<u8> {
    let mut out = Vec::with_capacity(buf.len());
    for line in buf.split(|b| *b == b'\n') {
        if line.is_empty() {
            out.push(b'\n');
            continue;
        }
        if line[0] == b'#' {
            out.extend_from_slice(line);
            out.push(b'\n');
            continue;
        }
        let s = match std::str::from_utf8(line) {
            Ok(s) => s,
            Err(_) => {
                out.extend_from_slice(line);
                out.push(b'\n');
                continue;
            }
        };
        let mut cols = s.splitn(3, '\t');
        let chrom = match cols.next() {
            Some(c) => c,
            None => continue,
        };
        let pos_1b: u64 = match cols.next().and_then(|p| p.parse().ok()) {
            Some(p) => p,
            None => {
                out.extend_from_slice(line);
                out.push(b'\n');
                continue;
            }
        };
        let pos_0b = pos_1b.saturating_sub(1);
        if filter.contains_position(chrom, pos_0b) {
            out.extend_from_slice(line);
            out.push(b'\n');
        }
    }
    if out.ends_with(b"\n\n") {
        out.pop();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixed_lengths() -> HashMap<String, u64> {
        let mut m = HashMap::new();
        m.insert("chr20".to_string(), 64_444_167);
        m.insert("chr21".to_string(), 46_709_983);
        m
    }

    #[test]
    fn empty_filter_is_unrestricted() {
        let f = RegionFilter::empty();
        assert!(f.is_unrestricted());
        assert!(!f.contains_position("chr20", 100));

        assert!(f.intersect_chunk_1based("chr20", 1, 100).is_empty());
    }

    #[test]
    fn bed_parse_basic() {
        let text = "chr20\t100\t200\nchr20\t250\t300\nchr21\t0\t50\n";
        let f = RegionFilter::from_bed_text(text).unwrap();
        assert!(!f.is_unrestricted());
        assert!(f.contains_position("chr20", 100));
        assert!(f.contains_position("chr20", 199));
        assert!(!f.contains_position("chr20", 200));
        assert!(f.contains_position("chr20", 250));
        assert!(!f.contains_position("chr20", 220));
        assert!(f.contains_position("chr21", 49));
        assert!(!f.contains_position("chr21", 50));
        assert!(!f.contains_position("chr22", 100));
    }

    #[test]
    fn bed_parse_skips_comments() {
        let text = "# comment\ntrack name=foo\nbrowser pos chr1:1-100\nchr1\t10\t20\n";
        let f = RegionFilter::from_bed_text(text).unwrap();
        assert_eq!(f.intervals_for("chr1").unwrap().len(), 1);
    }

    #[test]
    fn bed_parse_rejects_short_lines() {
        let res = RegionFilter::from_bed_text("chr1\t100\n");
        assert!(res.is_err());
    }

    #[test]
    fn bed_merges_overlapping() {

        let text = "chr1\t100\t200\nchr1\t150\t250\nchr1\t300\t400\n";
        let f = RegionFilter::from_bed_text(text).unwrap();
        let ivs = f.intervals_for("chr1").unwrap();
        assert_eq!(ivs.len(), 2);
        assert_eq!(ivs[0], Interval0Based::new(100, 250));
        assert_eq!(ivs[1], Interval0Based::new(300, 400));
    }

    #[test]
    fn bed_merges_touching() {

        let text = "chr1\t100\t200\nchr1\t200\t300\n";
        let f = RegionFilter::from_bed_text(text).unwrap();
        let ivs = f.intervals_for("chr1").unwrap();
        assert_eq!(ivs.len(), 1);
        assert_eq!(ivs[0], Interval0Based::new(100, 300));
    }

    #[test]
    fn region_strings_parse() {
        let specs = ["chr20:100-200", "chr20:300-400"];
        let lens = fixed_lengths();
        let f = RegionFilter::from_region_strings(&specs, &lens).unwrap();

        let ivs = f.intervals_for("chr20").unwrap();
        assert_eq!(ivs.len(), 2);
        assert_eq!(ivs[0], Interval0Based::new(99, 200));
        assert_eq!(ivs[1], Interval0Based::new(299, 400));
    }

    #[test]
    fn region_strings_bare_contig_uses_full_length() {
        let lens = fixed_lengths();
        let f = RegionFilter::from_region_strings(&["chr20"], &lens).unwrap();
        let ivs = f.intervals_for("chr20").unwrap();
        assert_eq!(ivs.len(), 1);
        assert_eq!(ivs[0], Interval0Based::new(0, 64_444_167));
    }

    #[test]
    fn region_strings_bare_contig_unknown_errors() {
        let lens = fixed_lengths();
        let res = RegionFilter::from_region_strings(&["chrZ"], &lens);
        assert!(res.is_err());
    }

    #[test]
    fn region_strings_rejects_zero_or_inverted() {
        let lens = fixed_lengths();
        assert!(RegionFilter::from_region_strings(&["chr20:0-100"], &lens).is_err());
        assert!(RegionFilter::from_region_strings(&["chr20:200-100"], &lens).is_err());
    }

    #[test]
    fn region_strings_strips_commas() {
        let lens = fixed_lengths();
        let f = RegionFilter::from_region_strings(&["chr20:1,000-2,000"], &lens).unwrap();
        let ivs = f.intervals_for("chr20").unwrap();
        assert_eq!(ivs[0], Interval0Based::new(999, 2000));
    }

    #[test]
    fn intersect_chunk_fully_inside() {

        let mut f = RegionFilter::empty();
        f.push_raw("chr20", 100, 200);
        f.finalize();
        let out = f.intersect_chunk_1based("chr20", 50, 150);
        assert_eq!(out, vec![(101, 150)]);
    }

    #[test]
    fn intersect_chunk_disjoint_returns_empty() {
        let mut f = RegionFilter::empty();
        f.push_raw("chr20", 100, 200);
        f.finalize();

        let out = f.intersect_chunk_1based("chr20", 300, 400);
        assert!(out.is_empty());
    }

    #[test]
    fn intersect_chunk_unknown_contig_returns_empty() {
        let mut f = RegionFilter::empty();
        f.push_raw("chr20", 100, 200);
        f.finalize();
        let out = f.intersect_chunk_1based("chrZ", 1, 1_000_000);
        assert!(out.is_empty());
    }

    #[test]
    fn intersect_chunk_multiple_sub_intervals() {

        let mut f = RegionFilter::empty();
        f.push_raw("chr20", 100, 200);
        f.push_raw("chr20", 300, 400);
        f.finalize();
        let out = f.intersect_chunk_1based("chr20", 1, 500);
        assert_eq!(out, vec![(101, 200), (301, 400)]);
    }

    #[test]
    fn union_combines_two_filters() {
        let mut a = RegionFilter::empty();
        a.push_raw("chr20", 100, 200);
        a.finalize();
        let mut b = RegionFilter::empty();
        b.push_raw("chr20", 150, 250);
        b.push_raw("chr21", 0, 50);
        b.finalize();
        let u = a.union(b);
        assert_eq!(
            u.intervals_for("chr20").unwrap(),
            &[Interval0Based::new(100, 250)]
        );
        assert_eq!(
            u.intervals_for("chr21").unwrap(),
            &[Interval0Based::new(0, 50)]
        );
    }

    #[test]
    fn covered_bp() {
        let mut f = RegionFilter::empty();
        f.push_raw("chr20", 100, 200);
        f.push_raw("chr20", 300, 400);
        f.finalize();
        assert_eq!(f.covered_bp("chr20"), 200);
        assert_eq!(f.covered_bp("chr21"), 0);
    }

    #[test]
    fn filter_vcf_preserves_headers_and_filters_rows() {
        let mut f = RegionFilter::empty();
        f.push_raw("chr20", 100, 200);
        f.finalize();
        let vcf = b"##fileformat=VCFv4.2\n\
                    ##source=clair3-connect\n\
                    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
                    chr20\t50\t.\tA\tG\t30\tPASS\t.\n\
                    chr20\t150\t.\tA\tC\t30\tPASS\t.\n\
                    chr20\t250\t.\tA\tT\t30\tPASS\t.\n\
                    chr21\t150\t.\tA\tT\t30\tPASS\t.\n";
        let out = filter_vcf_by_region(vcf, &f);
        let text = std::str::from_utf8(&out).unwrap();
        assert!(text.contains("##fileformat=VCFv4.2"));
        assert!(text.contains("##source=clair3-connect"));
        assert!(text.contains("#CHROM"));

        assert!(!text.contains("\t50\t."));

        assert!(text.contains("\t150\t."));

        assert!(!text.contains("\t250\t."));

        assert!(!text.contains("chr21\t150"));
    }

    #[test]
    fn filter_vcf_passes_through_malformed_rows() {

        let mut f = RegionFilter::empty();
        f.push_raw("chr20", 100, 200);
        f.finalize();
        let vcf = b"#CHROM\tPOS\nchr20\tNOT_A_NUMBER\textra\n";
        let out = filter_vcf_by_region(vcf, &f);
        let text = std::str::from_utf8(&out).unwrap();
        assert!(text.contains("NOT_A_NUMBER"));
    }

    #[test]
    fn contains_position_binary_search_correctness() {
        let mut f = RegionFilter::empty();

        for i in 0..1000 {
            let s = (i as u64) * 100;
            f.push_raw("chr1", s, s + 50);
        }
        f.finalize();

        assert!(f.contains_position("chr1", 25));

        assert!(!f.contains_position("chr1", 75));

        assert!(!f.contains_position("chr1", 99_975));

        assert!(f.contains_position("chr1", 99_925));
    }
}
