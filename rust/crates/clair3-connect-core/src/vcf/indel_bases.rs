#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InsertionLookupResult {

    Single(String),

    Multi(Vec<String>),
}

impl InsertionLookupResult {

    pub fn is_empty(&self) -> bool {
        match self {
            Self::Single(s) => s.is_empty(),
            Self::Multi(v) => v.is_empty(),
        }
    }

    pub fn best(&self) -> &str {
        match self {
            Self::Single(s) => s.as_str(),
            Self::Multi(v) => v.first().map(|s| s.as_str()).unwrap_or(""),
        }
    }
}

pub type DeletionLookupResult = InsertionLookupResult;

#[derive(Debug, Clone, Copy)]
pub struct InsertionOptions<'a> {

    pub propose_insertion_length: Option<usize>,

    pub minimum_insertion_length: usize,

    pub maximum_insertion_length: usize,

    pub insertion_bases_to_ignore: &'a str,

    pub return_multi: bool,
}

impl<'a> Default for InsertionOptions<'a> {
    fn default() -> Self {
        Self {
            propose_insertion_length: None,
            minimum_insertion_length: 1,
            maximum_insertion_length: 50,
            insertion_bases_to_ignore: "",
            return_multi: false,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct DeletionOptions<'a> {

    pub propose_deletion_length: Option<usize>,

    pub minimum_deletion_length: usize,

    pub maximum_deletion_length: usize,

    pub deletion_bases_to_ignore: &'a str,

    pub return_multi: bool,
}

impl<'a> Default for DeletionOptions<'a> {
    fn default() -> Self {
        Self {
            propose_deletion_length: None,
            minimum_deletion_length: 1,
            maximum_deletion_length: 50,
            deletion_bases_to_ignore: "",
            return_multi: false,
        }
    }
}

pub fn insertion_bases_using_alt_info_from(
    entries: &[(String, u32)],
    opts: InsertionOptions<'_>,
) -> InsertionLookupResult {
    if entries.is_empty() {
        return InsertionLookupResult::Single(String::new());
    }

    let propose_with_refbase: Option<usize> = opts.propose_insertion_length.map(|n| n + 1);

    let mut propose_matches: Vec<(&str, u32)> = Vec::new();
    let mut in_range: Vec<(&str, u32)> = Vec::new();

    for (raw_key, count) in entries {
        let bytes = raw_key.as_bytes();
        if bytes.is_empty() || bytes[0] != b'I' {
            continue;
        }

        let key = &raw_key[1..];
        if key == opts.insertion_bases_to_ignore {
            continue;
        }
        let key_len = key.len();

        if let Some(target) = propose_with_refbase {
            if key_len == target {
                propose_matches.push((key, *count));
            }
        }
        if key_len >= opts.minimum_insertion_length && key_len <= opts.maximum_insertion_length {
            in_range.push((key, *count));
        }
    }

    if propose_with_refbase.is_some() && !propose_matches.is_empty() {
        let best = argmax_first_seen(&propose_matches);
        return InsertionLookupResult::Single(best.to_string());
    }

    if opts.return_multi {
        let mut sorted = in_range.clone();

        sorted.sort_by_key(|(_, c)| *c);
        let topn: Vec<String> = sorted
            .iter()
            .rev()
            .take(2)
            .map(|(s, _)| (*s).to_string())
            .collect();
        return InsertionLookupResult::Multi(topn);
    }

    if in_range.is_empty() {
        return InsertionLookupResult::Single(String::new());
    }
    let best = argmax_first_seen(&in_range);
    InsertionLookupResult::Single(best.to_string())
}

pub fn deletion_bases_using_alt_info_from(
    entries: &[(String, u32)],
    opts: DeletionOptions<'_>,
) -> DeletionLookupResult {
    if entries.is_empty() {
        return DeletionLookupResult::Single(String::new());
    }

    let mut propose_matches: Vec<(&str, u32)> = Vec::new();
    let mut in_range: Vec<(&str, u32)> = Vec::new();

    for (raw_key, count) in entries {
        let bytes = raw_key.as_bytes();
        if bytes.is_empty() || bytes[0] != b'D' {
            continue;
        }
        let key = &raw_key[1..];
        if key == opts.deletion_bases_to_ignore {
            continue;
        }
        let key_len = key.len();

        if let Some(target) = opts.propose_deletion_length {
            if key_len == target {
                propose_matches.push((key, *count));
            }
        }
        if key_len >= opts.minimum_deletion_length && key_len <= opts.maximum_deletion_length {
            in_range.push((key, *count));
        }
    }

    if opts.propose_deletion_length.is_some() && !propose_matches.is_empty() {
        let best = argmax_first_seen(&propose_matches);
        return DeletionLookupResult::Single(best.to_string());
    }

    if opts.return_multi {
        let mut sorted = in_range.clone();
        sorted.sort_by_key(|(_, c)| *c);

        let top_two: Vec<&str> = sorted.iter().rev().take(2).map(|(s, _)| *s).collect();
        if top_two.len() < 2 {
            return DeletionLookupResult::Multi(Vec::new());
        }

        let (a, b) = (top_two[0], top_two[1]);
        let ordered = if a.len() > b.len() {
            vec![a.to_string(), b.to_string()]
        } else {
            vec![b.to_string(), a.to_string()]
        };
        return DeletionLookupResult::Multi(ordered);
    }

    if in_range.is_empty() {
        return DeletionLookupResult::Single(String::new());
    }
    let best = argmax_first_seen(&in_range);
    DeletionLookupResult::Single(best.to_string())
}

fn argmax_first_seen<'a>(items: &'a [(&'a str, u32)]) -> &'a str {
    debug_assert!(!items.is_empty());
    let mut best_idx = 0;
    let mut best_count = items[0].1;
    for (i, (_, c)) in items.iter().enumerate().skip(1) {
        if *c > best_count {
            best_idx = i;
            best_count = *c;
        }
    }
    items[best_idx].0
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::find_alt_base::AltInfoEntries;

    fn entries(pairs: &[(&str, u32)]) -> AltInfoEntries {
        pairs.iter().map(|(k, v)| (k.to_string(), *v)).collect()
    }

    #[test]
    fn ins_empty_returns_empty_string() {
        let r = insertion_bases_using_alt_info_from(&entries(&[]), InsertionOptions::default());
        assert_eq!(r, InsertionLookupResult::Single(String::new()));
        assert!(r.is_empty());
    }

    #[test]
    fn ins_argmax_picks_highest_support() {

        let r = insertion_bases_using_alt_info_from(
            &entries(&[("IAT", 6), ("IAA", 3), ("RA", 30), ("XG", 1)]),
            InsertionOptions::default(),
        );
        assert_eq!(r, InsertionLookupResult::Single("AT".into()));
    }

    #[test]
    fn ins_proposed_length_filters() {

        let r = insertion_bases_using_alt_info_from(
            &entries(&[("IAT", 30), ("IAAT", 5)]),
            InsertionOptions {
                propose_insertion_length: Some(2),
                ..Default::default()
            },
        );
        assert_eq!(r, InsertionLookupResult::Single("AAT".into()));
    }

    #[test]
    fn ins_proposed_length_no_match_falls_through_to_argmax() {

        let r = insertion_bases_using_alt_info_from(
            &entries(&[("IAT", 5), ("IAAT", 9)]),
            InsertionOptions {
                propose_insertion_length: Some(99),
                ..Default::default()
            },
        );
        assert_eq!(r, InsertionLookupResult::Single("AAT".into()));
    }

    #[test]
    fn ins_bases_to_ignore_skipped() {
        let r = insertion_bases_using_alt_info_from(
            &entries(&[("IAA", 9), ("IAT", 5)]),
            InsertionOptions {
                insertion_bases_to_ignore: "AA",
                ..Default::default()
            },
        );
        assert_eq!(r, InsertionLookupResult::Single("AT".into()));
    }

    #[test]
    fn ins_multi_returns_top_two_count_desc() {
        let r = insertion_bases_using_alt_info_from(
            &entries(&[("IAA", 3), ("IAT", 8), ("IAC", 5), ("IAG", 1)]),
            InsertionOptions {
                return_multi: true,
                ..Default::default()
            },
        );

        assert_eq!(
            r,
            InsertionLookupResult::Multi(vec!["AT".into(), "AC".into()])
        );
    }

    #[test]
    fn ins_skips_long_entries() {

        let long_ins = format!("I{}", "A".repeat(60));
        let r = insertion_bases_using_alt_info_from(
            &entries(&[(long_ins.as_str(), 100), ("IAA", 1)]),
            InsertionOptions::default(),
        );
        assert_eq!(r, InsertionLookupResult::Single("AA".into()));
    }

    #[test]
    fn ins_first_seen_wins_on_count_tie() {

        let r = insertion_bases_using_alt_info_from(
            &entries(&[("IAC", 5), ("IAT", 5), ("IAA", 5)]),
            InsertionOptions::default(),
        );
        assert_eq!(r, InsertionLookupResult::Single("AC".into()));
    }

    #[test]
    fn del_argmax() {

        let r = deletion_bases_using_alt_info_from(
            &entries(&[("DT", 5), ("DTT", 3)]),
            DeletionOptions::default(),
        );
        assert_eq!(r, DeletionLookupResult::Single("T".into()));
    }

    #[test]
    fn del_proposed_length_no_plus_one() {

        let r = deletion_bases_using_alt_info_from(
            &entries(&[("DT", 99), ("DTT", 1), ("DTTT", 50)]),
            DeletionOptions {
                propose_deletion_length: Some(2),
                ..Default::default()
            },
        );
        assert_eq!(r, DeletionLookupResult::Single("TT".into()));
    }

    #[test]
    fn del_multi_orders_longer_first() {

        let r = deletion_bases_using_alt_info_from(
            &entries(&[("DT", 8), ("DTT", 5), ("DTTT", 3)]),
            DeletionOptions {
                return_multi: true,
                ..Default::default()
            },
        );
        assert_eq!(
            r,
            DeletionLookupResult::Multi(vec!["TT".into(), "T".into()])
        );
    }

    #[test]
    fn del_multi_returns_empty_when_only_one_deletion() {
        let r = deletion_bases_using_alt_info_from(
            &entries(&[("DT", 8)]),
            DeletionOptions {
                return_multi: true,
                ..Default::default()
            },
        );

        assert_eq!(r, DeletionLookupResult::Multi(Vec::new()));
    }

    #[test]
    fn del_bases_to_ignore_skipped() {
        let r = deletion_bases_using_alt_info_from(
            &entries(&[("DT", 9), ("DTT", 5)]),
            DeletionOptions {
                deletion_bases_to_ignore: "T",
                ..Default::default()
            },
        );
        assert_eq!(r, DeletionLookupResult::Single("TT".into()));
    }

    #[test]
    fn ignores_non_indel_entries_in_both() {
        let e = entries(&[("XG", 31), ("RA", 2), ("IAT", 5), ("DT", 3)]);
        let r_ins = insertion_bases_using_alt_info_from(&e, InsertionOptions::default());
        let r_del = deletion_bases_using_alt_info_from(&e, DeletionOptions::default());
        assert_eq!(r_ins, InsertionLookupResult::Single("AT".into()));
        assert_eq!(r_del, DeletionLookupResult::Single("T".into()));
    }
}
