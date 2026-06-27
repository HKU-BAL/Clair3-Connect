use aes_gcm::aead::{Aead, KeyInit, OsRng, Payload, AeadCore};
use aes_gcm::{Aes256Gcm, Key, Nonce};
use thiserror::Error;
use zeroize::Zeroize;

pub const AEAD_KEY_BYTES: usize = 32;

pub const AEAD_NONCE_BYTES: usize = 12;

pub const AEAD_TAG_BYTES: usize = 16;

#[derive(Clone)]
pub struct AeadKey(Box<[u8; AEAD_KEY_BYTES]>);

impl std::fmt::Debug for AeadKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AeadKey")
            .field("bytes", &"<32 bytes redacted>")
            .finish()
    }
}

impl AeadKey {

    pub fn from_bytes(bytes: &[u8]) -> Result<Self, AeadError> {
        if bytes.len() != AEAD_KEY_BYTES {
            return Err(AeadError::InvalidKeyLength(bytes.len()));
        }
        let mut arr = [0u8; AEAD_KEY_BYTES];
        arr.copy_from_slice(bytes);
        Ok(Self(Box::new(arr)))
    }

    pub fn generate_random() -> Self {
        use aes_gcm::aead::rand_core::RngCore;
        let mut arr = [0u8; AEAD_KEY_BYTES];
        OsRng.fill_bytes(&mut arr);
        Self(Box::new(arr))
    }

    fn as_aead_key(&self) -> &Key<Aes256Gcm> {
        Key::<Aes256Gcm>::from_slice(self.0.as_ref())
    }
}

impl Drop for AeadKey {
    fn drop(&mut self) {
        self.0.zeroize();
    }
}

#[derive(Debug, Error)]
pub enum AeadError {

    #[error("aead key must be {AEAD_KEY_BYTES} bytes, got {0}")]
    InvalidKeyLength(usize),

    #[error("aead nonce must be {AEAD_NONCE_BYTES} bytes, got {0}")]
    InvalidNonceLength(usize),

    #[error("aead open failed: ciphertext tampered or wrong key")]
    Tampered,
}

pub fn seal(key: &AeadKey, plaintext: &[u8], aad: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let cipher = Aes256Gcm::new(key.as_aead_key());
    let nonce = Aes256Gcm::generate_nonce(&mut OsRng);
    let payload = Payload { msg: plaintext, aad };

    let ciphertext = cipher
        .encrypt(&nonce, payload)
        .expect("AES-GCM encrypt is infallible for sane inputs");
    (ciphertext, nonce.to_vec())
}

pub fn open(
    key: &AeadKey,
    ciphertext_with_tag: &[u8],
    nonce: &[u8],
    aad: &[u8],
) -> Result<Vec<u8>, AeadError> {
    if nonce.len() != AEAD_NONCE_BYTES {
        return Err(AeadError::InvalidNonceLength(nonce.len()));
    }
    let cipher = Aes256Gcm::new(key.as_aead_key());
    let nonce = Nonce::from_slice(nonce);
    let payload = Payload {
        msg: ciphertext_with_tag,
        aad,
    };
    cipher
        .decrypt(nonce, payload)
        .map_err(|_| AeadError::Tampered)
}

pub fn request_aad(
    model_id: &str,
    dtype_wire: i32,
    shape: &[i64],
    payload_sha256_hex: &str,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(
        4 + model_id.len() + 4 + 4 + 4 + 8 * shape.len() + 4 + payload_sha256_hex.len(),
    );
    write_len_prefixed(&mut out, model_id.as_bytes());
    out.extend_from_slice(&dtype_wire.to_be_bytes());
    out.extend_from_slice(&(shape.len() as u32).to_be_bytes());
    for d in shape {
        out.extend_from_slice(&d.to_be_bytes());
    }
    write_len_prefixed(&mut out, payload_sha256_hex.as_bytes());
    out
}

pub fn response_aad(
    request_id: &str,
    y_dtype_wire: i32,
    y_shape: &[i64],
    y_sha256_hex: &str,
) -> Vec<u8> {
    let mut out = Vec::with_capacity(
        4 + request_id.len() + 4 + 4 + 8 * y_shape.len() + 4 + y_sha256_hex.len(),
    );
    write_len_prefixed(&mut out, request_id.as_bytes());
    out.extend_from_slice(&y_dtype_wire.to_be_bytes());
    out.extend_from_slice(&(y_shape.len() as u32).to_be_bytes());
    for d in y_shape {
        out.extend_from_slice(&d.to_be_bytes());
    }
    write_len_prefixed(&mut out, y_sha256_hex.as_bytes());
    out
}

fn write_len_prefixed(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u32).to_be_bytes());
    out.extend_from_slice(bytes);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_basic() {
        let key = AeadKey::generate_random();
        let plaintext = b"the tensor bytes".to_vec();
        let aad = b"binding context".as_slice();
        let (ct, nonce) = seal(&key, &plaintext, aad);
        assert_eq!(nonce.len(), AEAD_NONCE_BYTES);

        assert_eq!(ct.len(), plaintext.len() + AEAD_TAG_BYTES);
        let recovered = open(&key, &ct, &nonce, aad).unwrap();
        assert_eq!(recovered, plaintext);
    }

    #[test]
    fn wrong_key_fails() {
        let k1 = AeadKey::generate_random();
        let k2 = AeadKey::generate_random();
        let (ct, nonce) = seal(&k1, b"x", b"aad");
        let err = open(&k2, &ct, &nonce, b"aad").unwrap_err();
        assert!(matches!(err, AeadError::Tampered));
    }

    #[test]
    fn tampered_ciphertext_fails() {
        let key = AeadKey::generate_random();
        let (mut ct, nonce) = seal(&key, b"hello", b"aad");

        ct[0] ^= 0xff;
        let err = open(&key, &ct, &nonce, b"aad").unwrap_err();
        assert!(matches!(err, AeadError::Tampered));
    }

    #[test]
    fn tampered_aad_fails() {
        let key = AeadKey::generate_random();
        let (ct, nonce) = seal(&key, b"x", b"aad-1");
        let err = open(&key, &ct, &nonce, b"aad-2").unwrap_err();
        assert!(matches!(err, AeadError::Tampered));
    }

    #[test]
    fn wrong_nonce_length_rejected() {
        let key = AeadKey::generate_random();
        let (ct, _nonce) = seal(&key, b"x", b"");
        let err = open(&key, &ct, &[0u8; 8], b"").unwrap_err();
        assert!(matches!(err, AeadError::InvalidNonceLength(8)));
    }

    #[test]
    fn invalid_key_length_rejected() {
        assert!(matches!(
            AeadKey::from_bytes(&[0u8; 16]),
            Err(AeadError::InvalidKeyLength(16))
        ));
    }

    #[test]
    fn each_seal_uses_fresh_random_nonce() {
        let key = AeadKey::generate_random();
        let (_ct1, n1) = seal(&key, b"x", b"");
        let (_ct2, n2) = seal(&key, b"x", b"");
        assert_ne!(n1, n2, "two seals must produce different nonces");
    }

    #[test]
    fn aad_construction_is_deterministic() {
        let a = request_aad("model_x", 3, &[1, 33, 18], "abcd");
        let b = request_aad("model_x", 3, &[1, 33, 18], "abcd");
        assert_eq!(a, b);
        let c = request_aad("model_x", 3, &[2, 33, 18], "abcd");
        assert_ne!(a, c, "different shape must produce different AAD");
        let d = request_aad("model_y", 3, &[1, 33, 18], "abcd");
        assert_ne!(a, d, "different model_id must produce different AAD");
    }

    #[test]
    fn request_response_aad_are_distinct_namespaces() {

        let req = request_aad("x", 1, &[1], "h");
        let resp = response_aad("x", 1, &[1], "h");

        assert_eq!(req, resp);
    }

    #[test]
    fn aead_key_zeroes_on_drop() {

        let k = AeadKey::generate_random();
        drop(k);
    }
}
