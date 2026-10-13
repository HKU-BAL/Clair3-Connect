use clair3_connect_transport::{aead, codec};

#[test]
fn compress_seal_open_decompress_roundtrip() {
    let mut tensor = vec![0u8; 32 * 89 * 33 * 8];
    for (i, b) in tensor.iter_mut().enumerate() {
        if i % 97 == 0 {
            *b = (i % 251) as u8;
        }
    }
    let key = aead::AeadKey::from_bytes(&[7u8; 32]).expect("key");
    let aad = b"model|dtype|shape|sha";

    let compressed = codec::compress_payload(&tensor).expect("compress");
    assert!(compressed.len() * 5 < tensor.len(), "tensor should compress several-fold");
    let (ciphertext, nonce) = aead::seal(&key, &compressed, aad);
    assert!(ciphertext.len() < tensor.len() / 5);

    let opened = aead::open(&key, &ciphertext, &nonce, aad).expect("open");
    let plain = codec::decompress_payload(&opened, tensor.len()).expect("decompress");
    assert_eq!(plain, tensor);

    assert!(codec::decompress_payload(&opened, tensor.len() / 2).is_err());
}
