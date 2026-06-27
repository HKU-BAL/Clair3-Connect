pub fn iupac_to_acgt(c: u8) -> u8 {
    match c {
        b'A' => b'A',
        b'C' => b'C',
        b'G' => b'G',
        b'T' => b'T',
        b'U' => b'T',
        b'R' => b'A',
        b'Y' => b'C',
        b'S' => b'C',
        b'W' => b'A',
        b'K' => b'G',
        b'M' => b'A',
        b'B' => b'C',
        b'D' => b'A',
        b'H' => b'A',
        b'V' => b'A',
        b'N' => b'A',

        b'a' => b'A',
        b'c' => b'C',
        b'g' => b'G',
        b't' | b'u' => b'T',
        _ => b'A',
    }
}

pub fn is_basic_base(c: u8) -> bool {
    matches!(c, b'A' | b'C' | b'G' | b'T' | b'U')
}

pub fn convert_iupac_to_n(s: &str) -> String {
    if s == "." {
        return s.to_string();
    }
    let mut out = String::with_capacity(s.len());
    let mut not_acgt = 0usize;
    for c in s.chars() {

        let upper = c.to_ascii_uppercase();
        if "ACGTN,.".contains(upper) {
            out.push(c);
        } else {
            not_acgt += 1;
            out.push('N');
        }
    }
    if not_acgt == 0 {
        s.to_string()
    } else {
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn acgt_passthrough() {
        for c in b"ACGT".iter().copied() {
            assert_eq!(iupac_to_acgt(c), c);
        }
    }

    #[test]
    fn iupac_collapse_table_matches_python() {

        let pairs = [
            (b'A', b'A'),
            (b'C', b'C'),
            (b'G', b'G'),
            (b'T', b'T'),
            (b'U', b'T'),
            (b'R', b'A'),
            (b'Y', b'C'),
            (b'S', b'C'),
            (b'W', b'A'),
            (b'K', b'G'),
            (b'M', b'A'),
            (b'B', b'C'),
            (b'D', b'A'),
            (b'H', b'A'),
            (b'V', b'A'),
            (b'N', b'A'),
        ];
        for (iupac, acgt) in pairs {
            assert_eq!(iupac_to_acgt(iupac), acgt, "iupac {}", iupac as char);
        }
    }

    #[test]
    fn basic_base_set() {
        assert!(is_basic_base(b'A'));
        assert!(is_basic_base(b'C'));
        assert!(is_basic_base(b'G'));
        assert!(is_basic_base(b'T'));
        assert!(is_basic_base(b'U'));
        assert!(!is_basic_base(b'N'));
        assert!(!is_basic_base(b'R'));
    }

    #[test]
    fn dot_passes_through() {
        assert_eq!(convert_iupac_to_n("."), ".");
    }

    #[test]
    fn all_acgtn_unchanged() {
        assert_eq!(convert_iupac_to_n("ACGTN"), "ACGTN");
        assert_eq!(convert_iupac_to_n("A,C"), "A,C");
    }

    #[test]
    fn non_acgt_becomes_n() {
        assert_eq!(convert_iupac_to_n("AYG"), "ANG");
        assert_eq!(convert_iupac_to_n("RYK"), "NNN");
        assert_eq!(convert_iupac_to_n("A.G"), "A.G");
    }

    #[test]
    fn early_exit_when_clean() {

        let s = "ACGTN";
        assert_eq!(convert_iupac_to_n(s), s);
    }
}
