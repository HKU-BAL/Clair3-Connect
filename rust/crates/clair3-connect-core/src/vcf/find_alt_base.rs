use thiserror::Error;

pub const MAX_DEPTH_GAP: u32 = 9;

pub type AltInfoEntries = Vec<(String, u32)>;

#[derive(Debug, Error)]
pub enum AltInfoError {

    #[error("alt_info payload missing '-' separator: {0:?}")]
    MissingSeparator(String),

    #[error("alt_info depth field not a u32: {0:?}")]
    InvalidDepth(String),

    #[error("alt_info count field not a u32 for tag {tag:?}: {raw:?}")]
    InvalidCount {

        tag: String,

        raw: String,
    },
}

pub fn find_alt_base(
    entries: &[(String, u32)],
    alternate_base: Option<char>,
) -> (Vec<char>, Option<char>) {

    let mut x_entries: Vec<(char, u32)> = entries
        .iter()
        .filter_map(|(k, v)| {
            let bytes = k.as_bytes();
            if bytes.len() >= 2 && bytes[0] == b'X' {
                Some((bytes[1] as char, *v))
            } else {
                None
            }
        })
        .collect();

    x_entries.sort_by_key(|e| std::cmp::Reverse(e.1));

    if x_entries.is_empty() {
        return (Vec::new(), None);
    }

    let top_count = x_entries[0].1;
    let alt_count_for_input = alternate_base.and_then(|ab| {
        x_entries
            .iter()
            .find(|(b, _)| *b == ab)
            .map(|(_, c)| *c)
    });

    let rewritten = match alt_count_for_input {
        None => Some(x_entries[0].0),
        Some(c) if top_count - c >= MAX_DEPTH_GAP => Some(x_entries[0].0),
        _ => alternate_base,
    };

    let sorted_alt_bases: Vec<char> = x_entries.iter().map(|(b, _)| *b).collect();
    (sorted_alt_bases, rewritten)
}

pub fn decode_alt_info(payload: &str) -> Result<(u32, AltInfoEntries), AltInfoError> {
    let (depth_str, rest) = payload
        .split_once('-')
        .ok_or_else(|| AltInfoError::MissingSeparator(payload.to_string()))?;
    let depth: u32 = depth_str
        .parse()
        .map_err(|_| AltInfoError::InvalidDepth(depth_str.to_string()))?;

    let toks: Vec<&str> = rest.split(' ').filter(|s| !s.is_empty()).collect();

    let pair_count = toks.len() / 2;
    let mut entries: AltInfoEntries = Vec::with_capacity(pair_count);
    for i in 0..pair_count {
        let tag = toks[i * 2];
        let raw_count = toks[i * 2 + 1];

        let count: u32 = raw_count.parse().map_err(|_| AltInfoError::InvalidCount {
            tag: tag.to_string(),
            raw: raw_count.to_string(),
        })?;
        entries.push((tag.to_string(), count));
    }
    Ok((depth, entries))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entries(pairs: &[(&str, u32)]) -> AltInfoEntries {
        pairs.iter().map(|(k, v)| (k.to_string(), *v)).collect()
    }

    #[test]
    fn empty_alt_info_returns_none() {
        let (sorted, alt) = find_alt_base(&entries(&[]), Some('G'));
        assert!(sorted.is_empty());
        assert!(alt.is_none());
    }

    #[test]
    fn rewrites_when_proposed_base_absent() {
        let (sorted, alt) = find_alt_base(&entries(&[("XC", 31), ("RA", 2)]), Some('G'));
        assert_eq!(sorted, vec!['C']);
        assert_eq!(alt, Some('C'));
    }

    #[test]
    fn keeps_proposed_within_gap() {
        let (sorted, alt) = find_alt_base(
            &entries(&[("XG", 31), ("XC", 24), ("RA", 2)]),
            Some('C'),
        );
        assert_eq!(sorted, vec!['G', 'C']);
        assert_eq!(alt, Some('C'));
    }

    #[test]
    fn rewrites_when_gap_exceeded() {
        let (sorted, alt) = find_alt_base(
            &entries(&[("XG", 31), ("XC", 20), ("RA", 2)]),
            Some('C'),
        );
        assert_eq!(sorted, vec!['G', 'C']);
        assert_eq!(alt, Some('G'));
    }

    #[test]
    fn ignores_non_x_entries() {
        let (sorted, alt) = find_alt_base(
            &entries(&[("XG", 31), ("RA", 2), ("DAT", 5), ("IAG", 1)]),
            Some('G'),
        );
        assert_eq!(sorted, vec!['G']);
        assert_eq!(alt, Some('G'));
    }

    #[test]
    fn no_proposed_base_returns_top() {
        let (sorted, alt) = find_alt_base(&entries(&[("XG", 31), ("XC", 20)]), None);
        assert_eq!(sorted, vec!['G', 'C']);
        assert_eq!(alt, Some('G'));
    }

    #[test]
    fn ties_resolved_by_input_order() {

        let (sorted, alt) = find_alt_base(
            &entries(&[("XC", 31), ("XG", 31), ("RA", 2)]),
            None,
        );
        assert_eq!(sorted, vec!['C', 'G'], "ties must preserve input order");
        assert_eq!(alt, Some('C'));
    }

    #[test]
    fn ties_other_input_order() {
        let (sorted, _alt) = find_alt_base(
            &entries(&[("XG", 31), ("XC", 31)]),
            None,
        );
        assert_eq!(sorted, vec!['G', 'C']);
    }

    #[test]
    fn decode_basic() {
        let (depth, e) = decode_alt_info("33-XG 31 RA 2").unwrap();
        assert_eq!(depth, 33);
        assert_eq!(e, vec![("XG".into(), 31), ("RA".into(), 2)]);
    }

    #[test]
    fn decode_indels() {
        let (depth, e) = decode_alt_info("28-DT 5 DTT 4 ICT 3 RC 4").unwrap();
        assert_eq!(depth, 28);
        assert_eq!(
            e,
            vec![
                ("DT".into(), 5),
                ("DTT".into(), 4),
                ("ICT".into(), 3),
                ("RC".into(), 4),
            ]
        );
    }

    #[test]
    fn decode_empty_payload() {
        let (depth, e) = decode_alt_info("0-").unwrap();
        assert_eq!(depth, 0);
        assert!(e.is_empty());
    }

    #[test]
    fn decode_rejects_missing_separator() {
        assert!(matches!(
            decode_alt_info("33"),
            Err(AltInfoError::MissingSeparator(_))
        ));
    }

    #[test]
    fn decode_rejects_non_numeric_depth() {
        assert!(matches!(
            decode_alt_info("abc-XG 1"),
            Err(AltInfoError::InvalidDepth(_))
        ));
    }

    #[test]
    fn decode_rejects_non_numeric_count() {
        assert!(matches!(
            decode_alt_info("33-XG NaN"),
            Err(AltInfoError::InvalidCount { .. })
        ));
    }

    #[test]
    fn decode_silently_drops_dangling_tag() {

        let (depth, e) = decode_alt_info("33-XG 31 XC").unwrap();
        assert_eq!(depth, 33);
        assert_eq!(e, vec![("XG".into(), 31)]);
    }
}
