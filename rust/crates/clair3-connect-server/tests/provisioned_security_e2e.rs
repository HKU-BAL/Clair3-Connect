use clair3_connect_server::provisioning::{generate_aead_key, generate_pki};
use clair3_connect_transport::aead::AeadKey;
use clair3_connect_transport::predict_client::{make_request_i8, PredictClient, RetryPolicy};
use clair3_connect_transport::proto::clair3_connect_inference_server::{
    Clair3ConnectInference, Clair3ConnectInferenceServer,
};
use clair3_connect_transport::proto::{
    HealthRequest, HealthResponse, PredictRequest, PredictResponse, PredictStatus,
};
use clair3_connect_transport::tls::{ClientTlsOptions, ServerTlsOptions};
use std::path::PathBuf;
use std::time::Duration;
use tokio::net::TcpListener;
use tonic::transport::Server;
use tonic::{Request, Response, Status};

#[derive(Default, Clone)]
struct EchoService;

#[tonic::async_trait]
impl Clair3ConnectInference for EchoService {
    async fn predict(
        &self,
        request: Request<PredictRequest>,
    ) -> Result<Response<PredictResponse>, Status> {
        let req = request.into_inner();
        Ok(Response::new(PredictResponse {
            request_id: req.request_id,
            status: PredictStatus::Ok as i32,
            ..Default::default()
        }))
    }

    async fn health(
        &self,
        _r: Request<HealthRequest>,
    ) -> Result<Response<HealthResponse>, Status> {
        Ok(Response::new(HealthResponse::default()))
    }
}

fn tmp_dir(label: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!(
        "clair3-connect-prov-e2e-{}-{}-{}",
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

async fn spawn_mtls_echo(server_tls: ServerTlsOptions) -> (String, tokio::task::JoinHandle<()>) {
    let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    let incoming = tokio_stream::wrappers::TcpListenerStream::new(listener);
    let handle = tokio::spawn(async move {
        Server::builder()
            .tls_config(server_tls.build())
            .unwrap()
            .add_service(Clair3ConnectInferenceServer::new(EchoService))
            .serve_with_incoming(incoming)
            .await
            .ok();
    });
    tokio::task::yield_now().await;
    (format!("https://localhost:{}", addr.port()), handle)
}

#[tokio::test]
async fn generated_pki_drives_mtls_roundtrip() {
    let dir = tmp_dir("ok");
    let pki = generate_pki(&dir, &["localhost".to_string(), "127.0.0.1".to_string()]).unwrap();

    let server_tls = ServerTlsOptions::from_paths(
        &pki.server_cert,
        &pki.server_key,
        Some(pki.ca_cert.as_path()),
    )
    .expect("server TLS from generated files");
    assert!(server_tls.is_mutual(), "client-CA present ⇒ mutual");
    let (endpoint, server) = spawn_mtls_echo(server_tls).await;

    let client_tls = ClientTlsOptions::from_paths(
        &pki.ca_cert,
        Some(pki.client_cert.as_path()),
        Some(pki.client_key.as_path()),
        "localhost",
    )
    .expect("client TLS from generated files");
    let policy = RetryPolicy {
        max_retries: 3,
        initial_delay: Duration::from_millis(50),
        max_delay: Duration::from_millis(200),
        per_attempt_timeout: None,
    };
    let mut client = PredictClient::connect_tls_with_retry(endpoint, &client_tls, policy)
        .await
        .expect("connect with generated client identity");

    let req = make_request_i8("prov-mtls", "echo-model", vec![1, 1], &[1]);
    let resp = client.predict(req).await.expect("predict over generated mTLS");
    assert_eq!(resp.request_id, "prov-mtls");
    assert_eq!(resp.status, PredictStatus::Ok as i32);

    server.abort();
    std::fs::remove_dir_all(&dir).ok();
}

#[tokio::test]
async fn generated_pki_rejects_client_without_identity() {
    let dir = tmp_dir("noident");
    let pki = generate_pki(&dir, &["localhost".to_string()]).unwrap();

    let server_tls = ServerTlsOptions::from_paths(
        &pki.server_cert,
        &pki.server_key,
        Some(pki.ca_cert.as_path()),
    )
    .unwrap();
    let (endpoint, server) = spawn_mtls_echo(server_tls).await;

    let client_tls = ClientTlsOptions::from_paths(&pki.ca_cert, None, None, "localhost").unwrap();
    let policy = RetryPolicy {
        max_retries: 1,
        initial_delay: Duration::from_millis(50),
        max_delay: Duration::from_millis(100),
        per_attempt_timeout: Some(Duration::from_secs(2)),
    };
    let outcome = match PredictClient::connect_tls_with_retry(endpoint, &client_tls, policy).await {
        Err(_) => true,
        Ok(mut c) => {
            let req = make_request_i8("noident", "echo-model", vec![1, 1], &[1]);
            c.predict(req).await.is_err()
        }
    };
    assert!(outcome, "client without identity must be rejected under mTLS");

    server.abort();
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn generated_aead_key_is_wire_loadable() {
    let dir = tmp_dir("aead");
    let key_path = dir.join("aead.key");
    generate_aead_key(&key_path).unwrap();

    let bytes = std::fs::read(&key_path).unwrap();
    let _key = AeadKey::from_bytes(&bytes).expect("generated AEAD key loads on the wire path");
    assert_eq!(bytes.len(), clair3_connect_transport::aead::AEAD_KEY_BYTES);
    std::fs::remove_dir_all(&dir).ok();
}
