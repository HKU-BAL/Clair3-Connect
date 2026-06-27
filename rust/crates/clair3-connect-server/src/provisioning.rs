use anyhow::{Context, Result};
use std::io::Read;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};

fn lock_down(path: &Path) -> Result<()> {
    let mut perms = std::fs::metadata(path)
        .with_context(|| format!("stat {}", path.display()))?
        .permissions();
    perms.set_mode(0o600);
    std::fs::set_permissions(path, perms)
        .with_context(|| format!("chmod 600 {}", path.display()))?;
    Ok(())
}

pub fn generate_aead_key(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("mkdir {}", parent.display()))?;
        }
    }
    let mut key = [0u8; clair3_connect_transport::aead::AEAD_KEY_BYTES];
    std::fs::File::open("/dev/urandom")
        .context("open /dev/urandom")?
        .read_exact(&mut key)
        .context("read 32 random bytes")?;
    std::fs::write(path, key).with_context(|| format!("write {}", path.display()))?;
    lock_down(path)?;
    Ok(())
}

#[derive(Debug, Clone)]
pub struct GeneratedPki {

    pub ca_cert: PathBuf,

    pub ca_key: PathBuf,

    pub server_cert: PathBuf,

    pub server_key: PathBuf,

    pub client_cert: PathBuf,

    pub client_key: PathBuf,
}

pub fn generate_pki(dir: &Path, sans: &[String]) -> Result<GeneratedPki> {
    use rcgen::{
        BasicConstraints, CertificateParams, DnType, IsCa, KeyPair, KeyUsagePurpose,
    };

    std::fs::create_dir_all(dir).with_context(|| format!("mkdir {}", dir.display()))?;

    let ca_key = KeyPair::generate().context("generate CA key")?;
    let mut ca_params =
        CertificateParams::new(Vec::<String>::new()).context("CA params")?;
    ca_params
        .distinguished_name
        .push(DnType::CommonName, "clair3-connect-local-CA");
    ca_params.is_ca = IsCa::Ca(BasicConstraints::Unconstrained);
    ca_params.key_usages = vec![
        KeyUsagePurpose::KeyCertSign,
        KeyUsagePurpose::DigitalSignature,
        KeyUsagePurpose::CrlSign,
    ];
    let ca_cert = ca_params.self_signed(&ca_key).context("self-sign CA")?;

    let server_key = KeyPair::generate().context("generate server key")?;
    let mut server_params =
        CertificateParams::new(sans.to_vec()).context("server params")?;
    server_params
        .distinguished_name
        .push(DnType::CommonName, "clair3-connect-server");
    let server_cert = server_params
        .signed_by(&server_key, &ca_cert, &ca_key)
        .context("sign server cert")?;

    let client_key = KeyPair::generate().context("generate client key")?;
    let mut client_params =
        CertificateParams::new(Vec::<String>::new()).context("client params")?;
    client_params
        .distinguished_name
        .push(DnType::CommonName, "clair3-connect-client");
    let client_cert = client_params
        .signed_by(&client_key, &ca_cert, &ca_key)
        .context("sign client cert")?;

    let out = GeneratedPki {
        ca_cert: dir.join("ca.crt"),
        ca_key: dir.join("ca.key"),
        server_cert: dir.join("server.crt"),
        server_key: dir.join("server.key"),
        client_cert: dir.join("client.crt"),
        client_key: dir.join("client.key"),
    };

    std::fs::write(&out.ca_cert, ca_cert.pem())
        .with_context(|| format!("write {}", out.ca_cert.display()))?;
    std::fs::write(&out.server_cert, server_cert.pem())
        .with_context(|| format!("write {}", out.server_cert.display()))?;
    std::fs::write(&out.client_cert, client_cert.pem())
        .with_context(|| format!("write {}", out.client_cert.display()))?;

    for (path, pem) in [
        (&out.ca_key, ca_key.serialize_pem()),
        (&out.server_key, server_key.serialize_pem()),
        (&out.client_key, client_key.serialize_pem()),
    ] {
        std::fs::write(path, pem).with_context(|| format!("write {}", path.display()))?;
        lock_down(path)?;
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_dir(label: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!(
            "clair3-connect-prov-{}-{}-{}",
            std::process::id(),
            label,
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn aead_key_is_32_bytes_and_locked() {
        let dir = tmp_dir("aead");
        let key = dir.join("aead.key");
        generate_aead_key(&key).unwrap();
        let bytes = std::fs::read(&key).unwrap();
        assert_eq!(bytes.len(), 32);

        clair3_connect_transport::aead::AeadKey::from_bytes(&bytes).unwrap();
        let mode = std::fs::metadata(&key).unwrap().permissions().mode();
        assert_eq!(mode & 0o777, 0o600);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn two_aead_keys_differ() {
        let dir = tmp_dir("aead2");
        let a = dir.join("a.key");
        let b = dir.join("b.key");
        generate_aead_key(&a).unwrap();
        generate_aead_key(&b).unwrap();
        assert_ne!(std::fs::read(&a).unwrap(), std::fs::read(&b).unwrap());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn pki_writes_six_pems_and_loads_as_tls() {
        let dir = tmp_dir("pki");
        let pki = generate_pki(&dir, &["localhost".to_string(), "127.0.0.1".to_string()])
            .unwrap();
        for p in [
            &pki.ca_cert,
            &pki.ca_key,
            &pki.server_cert,
            &pki.server_key,
            &pki.client_cert,
            &pki.client_key,
        ] {
            assert!(p.exists(), "missing {}", p.display());
            assert!(std::fs::read_to_string(p).unwrap().contains("-----BEGIN"));
        }

        for k in [&pki.ca_key, &pki.server_key, &pki.client_key] {
            assert_eq!(
                std::fs::metadata(k).unwrap().permissions().mode() & 0o777,
                0o600
            );
        }

        let opts = clair3_connect_transport::tls::ServerTlsOptions::from_paths(
            &pki.server_cert,
            &pki.server_key,
            Some(pki.ca_cert.as_path()),
        )
        .unwrap();
        assert!(opts.is_mutual());
        std::fs::remove_dir_all(&dir).ok();
    }
}
