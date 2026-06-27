use clair3_connect_transport::aead::{seal, AeadKey};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() != 5 {
        eprintln!("usage: aead_seal <key32> <plaintext_in> <ct1_out> <ct2_out>");
        std::process::exit(2);
    }
    let key_bytes = std::fs::read(&a[1]).expect("read key");
    let key = AeadKey::from_bytes(&key_bytes).expect("key must be 32 bytes");
    let pt = std::fs::read(&a[2]).expect("read plaintext");

    let aad = b"clair3-connect:full_alignment_ont_hac_v500:int8:[64,89,33,8]".as_slice();
    let (ct1, _n1) = seal(&key, &pt, aad);
    let (ct2, _n2) = seal(&key, &pt, aad);
    std::fs::write(&a[3], &ct1).expect("write ct1");
    std::fs::write(&a[4], &ct2).expect("write ct2");
    eprintln!(
        "sealed {} plaintext bytes -> {} ciphertext bytes (AES-256-GCM, +16B tag), x2 fresh nonces",
        pt.len(),
        ct1.len()
    );
}
