use std::fs;
use std::path::Path;

use tonic::transport::{Certificate, ClientTlsConfig, Identity, ServerTlsConfig};

#[derive(Debug, Clone)]
pub struct ServerTlsOptions {

    pub cert_pem: Vec<u8>,

    pub key_pem: Vec<u8>,

    pub client_ca_pem: Option<Vec<u8>>,
}

#[derive(Debug, Clone)]
pub struct ClientTlsOptions {

    pub server_ca_pem: Vec<u8>,

    pub client_identity_pem: Option<(Vec<u8>, Vec<u8>)>,

    pub domain: String,
}

impl ServerTlsOptions {

    pub fn from_paths(
        cert: &Path,
        key: &Path,
        client_ca: Option<&Path>,
    ) -> std::io::Result<Self> {
        let cert_pem = fs::read(cert)?;
        let key_pem = fs::read(key)?;
        let client_ca_pem = match client_ca {
            Some(p) => Some(fs::read(p)?),
            None => None,
        };
        Ok(Self {
            cert_pem,
            key_pem,
            client_ca_pem,
        })
    }

    pub fn build(&self) -> ServerTlsConfig {
        let identity = Identity::from_pem(&self.cert_pem, &self.key_pem);
        let mut cfg = ServerTlsConfig::new().identity(identity);
        if let Some(ca) = &self.client_ca_pem {
            cfg = cfg.client_ca_root(Certificate::from_pem(ca));
        }
        cfg
    }

    pub fn is_mutual(&self) -> bool {
        self.client_ca_pem.is_some()
    }
}

impl ClientTlsOptions {

    pub fn from_paths(
        server_ca: &Path,
        client_cert: Option<&Path>,
        client_key: Option<&Path>,
        domain: &str,
    ) -> std::io::Result<Self> {
        let server_ca_pem = fs::read(server_ca)?;
        let client_identity_pem = match (client_cert, client_key) {
            (Some(c), Some(k)) => Some((fs::read(c)?, fs::read(k)?)),
            (None, None) => None,
            _ => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "client_cert and client_key must be specified together",
                ))
            }
        };
        Ok(Self {
            server_ca_pem,
            client_identity_pem,
            domain: domain.to_string(),
        })
    }

    pub fn build(&self) -> ClientTlsConfig {
        let mut cfg = ClientTlsConfig::new()
            .ca_certificate(Certificate::from_pem(&self.server_ca_pem))
            .domain_name(&self.domain);
        if let Some((cert, key)) = &self.client_identity_pem {
            cfg = cfg.identity(Identity::from_pem(cert, key));
        }
        cfg
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_mutual_flag() {
        let opts = ServerTlsOptions {
            cert_pem: b"dummy".to_vec(),
            key_pem: b"dummy".to_vec(),
            client_ca_pem: None,
        };
        assert!(!opts.is_mutual());

        let opts = ServerTlsOptions {
            cert_pem: b"dummy".to_vec(),
            key_pem: b"dummy".to_vec(),
            client_ca_pem: Some(b"dummy".to_vec()),
        };
        assert!(opts.is_mutual());
    }

    #[test]
    fn from_paths_missing_file_errors() {
        let res = ServerTlsOptions::from_paths(
            Path::new("/nonexistent/cert.pem"),
            Path::new("/nonexistent/key.pem"),
            None,
        );
        assert!(res.is_err());
    }

    #[test]
    fn client_from_paths_partial_identity_errors() {
        let tmp = std::env::temp_dir();
        let ca = tmp.join("clair3_connect_test_ca.pem");
        std::fs::write(&ca, b"dummy ca").unwrap();
        let res = ClientTlsOptions::from_paths(&ca, Some(Path::new("/nope.pem")), None, "x");
        assert!(res.is_err());
        let res = ClientTlsOptions::from_paths(&ca, None, Some(Path::new("/nope.pem")), "x");
        assert!(res.is_err());
        std::fs::remove_file(&ca).ok();
    }
}
