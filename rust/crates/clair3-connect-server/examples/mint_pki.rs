use clair3_connect_server::provisioning::{generate_aead_key, generate_pki};
use std::path::Path;

fn main() {
    let dir = match std::env::args().nth(1) {
        Some(d) => d,
        None => {
            eprintln!("usage: mint_pki <dir>");
            std::process::exit(2);
        }
    };
    let dir = Path::new(&dir);
    let pki = generate_pki(dir, &["localhost".to_string(), "127.0.0.1".to_string()])
        .expect("generate_pki");
    let aead = dir.join("aead.key");
    generate_aead_key(&aead).expect("generate_aead_key");

    println!("ca_cert={}", pki.ca_cert.display());
    println!("ca_key={}", pki.ca_key.display());
    println!("server_cert={}", pki.server_cert.display());
    println!("server_key={}", pki.server_key.display());
    println!("client_cert={}", pki.client_cert.display());
    println!("client_key={}", pki.client_key.display());
    println!("aead_key={}", aead.display());
}
