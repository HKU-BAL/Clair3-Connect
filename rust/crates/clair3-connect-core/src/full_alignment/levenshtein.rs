pub fn levenshtein(a: &[u8], b: &[u8]) -> usize {

    if std::ptr::eq(a.as_ptr(), b.as_ptr()) && a.len() == b.len() {
        return 0;
    }
    let length = a.len();
    let b_length = b.len();
    if length == 0 {
        return b_length;
    }
    if b_length == 0 {
        return length;
    }

    let mut cache: Vec<usize> = (1..=length).collect();
    let mut result: usize = 0;
    for (b_index, &code) in b.iter().enumerate() {
        let mut distance = b_index;
        result = b_index;
        for (index, &a_byte) in a.iter().enumerate() {
            let b_distance = if code == a_byte { distance } else { distance + 1 };
            distance = cache[index];

            let new = if distance > result {
                if b_distance > result {
                    result + 1
                } else {
                    b_distance
                }
            } else if b_distance > distance {
                distance + 1
            } else {
                b_distance
            };
            cache[index] = new;
            result = new;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_inputs() {
        assert_eq!(levenshtein(b"", b""), 0);
        assert_eq!(levenshtein(b"", b"abc"), 3);
        assert_eq!(levenshtein(b"abc", b""), 3);
    }

    #[test]
    fn identical_strings_distance_zero() {
        assert_eq!(levenshtein(b"hello", b"hello"), 0);
        assert_eq!(levenshtein(b"AGCTACGT", b"AGCTACGT"), 0);
    }

    #[test]
    fn single_sub() {
        assert_eq!(levenshtein(b"AGCT", b"AGGT"), 1);
        assert_eq!(levenshtein(b"AAAA", b"ACAA"), 1);
    }

    #[test]
    fn single_insertion_or_deletion() {
        assert_eq!(levenshtein(b"AGCT", b"AGT"), 1);
        assert_eq!(levenshtein(b"AGT", b"AGCT"), 1);
    }

    #[test]
    fn classic_kitten_sitting() {

        assert_eq!(levenshtein(b"kitten", b"sitting"), 3);
    }

    #[test]
    fn snp_realign_example() {

        let query = b"AAAACATGCATGCATGCATG";
        let r_ref = b"AAAACATGCATGCATGCATG";
        let r_alt = b"AAAACATGCATCCATGCATG";
        assert_eq!(levenshtein(query, r_ref), 0);
        assert_eq!(levenshtein(query, r_alt), 1);
    }

    #[test]
    fn long_strings_dna() {
        let a = b"ACGTACGTACGTACGTACGTACGTACGT";
        let b = b"ACGTACGTACGTAGGTACGTACGTACGT";
        assert_eq!(levenshtein(a, b), 1);
    }
}
