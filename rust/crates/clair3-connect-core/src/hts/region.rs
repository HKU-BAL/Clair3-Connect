use anyhow::{anyhow, Context, Result};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {

    pub contig: String,

    pub start: u64,

    pub end: u64,
}

impl Region {

    pub fn parse(s: &str) -> Result<Self> {
        let (contig, rng) = s
            .split_once(':')
            .ok_or_else(|| anyhow!("region {s:?} missing ':'"))?;
        let (start_s, end_s) = rng
            .split_once('-')
            .ok_or_else(|| anyhow!("region {s:?} missing '-'"))?;
        let start_1: u64 = start_s
            .parse()
            .with_context(|| format!("parse start from {s:?}"))?;
        let end_1: u64 = end_s
            .parse()
            .with_context(|| format!("parse end from {s:?}"))?;
        if end_1 < start_1 {
            return Err(anyhow!("region {s:?} end < start"));
        }

        Ok(Self {
            contig: contig.to_string(),
            start: start_1.saturating_sub(1),
            end: end_1,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn parse_canonical() {
        let r = Region::parse("chr20:100000-300000").unwrap();
        assert_eq!(r.contig, "chr20");
        assert_eq!(r.start, 99_999);
        assert_eq!(r.end, 300_000);
    }
    #[test]
    fn parse_with_chr_prefix() {
        let r = Region::parse("chrM:1-100").unwrap();
        assert_eq!(r.contig, "chrM");
        assert_eq!(r.start, 0);
        assert_eq!(r.end, 100);
    }
    #[test]
    fn parse_missing_colon() {
        assert!(Region::parse("chr20-100000-300000").is_err());
    }
    #[test]
    fn parse_invalid_range() {
        assert!(Region::parse("chr20:300000-100000").is_err());
    }
}
