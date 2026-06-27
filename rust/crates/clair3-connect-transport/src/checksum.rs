use sha2::{Digest, Sha256};

pub fn sha256_hex(payload: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(payload);
    hex::encode(hasher.finalize())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_payload_canonical_digest() {
        assert_eq!(
            sha256_hex(&[]),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
    }

    #[test]
    fn deterministic_per_input() {
        let p = b"clair3-connect-test-payload";
        assert_eq!(sha256_hex(p), sha256_hex(p));
        assert_ne!(sha256_hex(p), sha256_hex(b"different"));
    }
}
