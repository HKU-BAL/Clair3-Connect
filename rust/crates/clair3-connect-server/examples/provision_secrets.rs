use clair3_connect_server::provisioning::{generate_aead_key, generate_pki};
use std::path::PathBuf;

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let dir = PathBuf::from(
        args.next()
            .expect("usage: provision_secrets <out_dir> [SAN ...]"),
    );
    let sans: Vec<String> = {
        let v: Vec<String> = args.collect();
        if v.is_empty() {
            vec!["localhost".to_string(), "127.0.0.1".to_string()]
        } else {
            v
        }
    };

    let pki = generate_pki(&dir, &sans)?;
    let aead = dir.join("aead.key");
    generate_aead_key(&aead)?;

    println!("Wrote PKI + AEAD key to {}", dir.display());
    println!(
        "server: --tls-cert {} --tls-key {} --tls-client-ca {} --aead-key {}",
        pki.server_cert.display(),
        pki.server_key.display(),
        pki.ca_cert.display(),
        aead.display()
    );
    println!(
        "client: --tls-server-ca {} --tls-client-cert {} --tls-client-key {} --aead-key {} --tls-domain {}",
        pki.ca_cert.display(),
        pki.client_cert.display(),
        pki.client_key.display(),
        aead.display(),
        sans[0]
    );
    Ok(())
}
