use clair3_connect_transport::predict_client::{
    make_request_i8, PredictClient, PredictError, RetryPolicy,
};
use clair3_connect_transport::proto::clair3_connect_inference_server::{
    Clair3ConnectInference, Clair3ConnectInferenceServer,
};
use clair3_connect_transport::proto::{
    HealthRequest, HealthResponse, PredictRequest, PredictResponse, PredictStatus,
};
use clair3_connect_transport::tls::{ClientTlsOptions, ServerTlsOptions};
use std::time::Duration;
use tokio::net::TcpListener;
use tonic::transport::Server;
use tonic::{Request, Response, Status};

struct TestCerts {
    ca_pem: Vec<u8>,

    other_ca_pem: Vec<u8>,
    server_cert_pem: Vec<u8>,
    server_key_pem: Vec<u8>,
    client_cert_pem: Vec<u8>,
    client_key_pem: Vec<u8>,
}

fn gen_ca_and_identities() -> TestCerts {
    use rcgen::{
        BasicConstraints, CertificateParams, DnType, IsCa, KeyPair, KeyUsagePurpose,
    };

    let ca_key = KeyPair::generate().unwrap();
    let mut ca_params = CertificateParams::new(Vec::<String>::new()).unwrap();
    ca_params
        .distinguished_name
        .push(DnType::CommonName, "clair3-connect-test-CA");
    ca_params.is_ca = IsCa::Ca(BasicConstraints::Unconstrained);
    ca_params.key_usages = vec![
        KeyUsagePurpose::KeyCertSign,
        KeyUsagePurpose::DigitalSignature,
        KeyUsagePurpose::CrlSign,
    ];
    let ca_cert = ca_params.self_signed(&ca_key).unwrap();

    let server_key = KeyPair::generate().unwrap();
    let mut server_params =
        CertificateParams::new(vec!["localhost".to_string()]).unwrap();
    server_params
        .distinguished_name
        .push(DnType::CommonName, "clair3-connect-test-server");
    let server_cert = server_params
        .signed_by(&server_key, &ca_cert, &ca_key)
        .unwrap();

    let client_key = KeyPair::generate().unwrap();
    let mut client_params = CertificateParams::new(Vec::<String>::new()).unwrap();
    client_params
        .distinguished_name
        .push(DnType::CommonName, "clair3-connect-test-client");
    let client_cert = client_params
        .signed_by(&client_key, &ca_cert, &ca_key)
        .unwrap();

    let other_ca_key = KeyPair::generate().unwrap();
    let mut other_ca_params = CertificateParams::new(Vec::<String>::new()).unwrap();
    other_ca_params
        .distinguished_name
        .push(DnType::CommonName, "unrelated-CA");
    other_ca_params.is_ca = IsCa::Ca(BasicConstraints::Unconstrained);
    let other_ca_cert = other_ca_params.self_signed(&other_ca_key).unwrap();

    TestCerts {
        ca_pem: ca_cert.pem().into_bytes(),
        other_ca_pem: other_ca_cert.pem().into_bytes(),
        server_cert_pem: server_cert.pem().into_bytes(),
        server_key_pem: server_key.serialize_pem().into_bytes(),
        client_cert_pem: client_cert.pem().into_bytes(),
        client_key_pem: client_key.serialize_pem().into_bytes(),
    }
}

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

async fn spawn_mtls_server(
    server_tls: ServerTlsOptions,
) -> (String, tokio::task::JoinHandle<()>) {
    let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    let incoming = tokio_stream::wrappers::TcpListenerStream::new(listener);
    let svc = EchoService;
    let handle = tokio::spawn(async move {
        Server::builder()
            .tls_config(server_tls.build())
            .unwrap()
            .add_service(
                Clair3ConnectInferenceServer::new(svc)
                    .accept_compressed(tonic::codec::CompressionEncoding::Zstd),
            )
            .serve_with_incoming(incoming)
            .await
            .ok();
    });

    tokio::task::yield_now().await;
    (format!("https://localhost:{}", addr.port()), handle)
}

#[tokio::test]
async fn mtls_predict_roundtrip() {
    let certs = gen_ca_and_identities();
    let server_tls = ServerTlsOptions {
        cert_pem: certs.server_cert_pem.clone(),
        key_pem: certs.server_key_pem.clone(),
        client_ca_pem: Some(certs.ca_pem.clone()),
    };
    let (endpoint, server) = spawn_mtls_server(server_tls).await;

    let client_tls = ClientTlsOptions {
        server_ca_pem: certs.ca_pem.clone(),
        client_identity_pem: Some((
            certs.client_cert_pem.clone(),
            certs.client_key_pem.clone(),
        )),
        domain: "localhost".to_string(),
    };
    let policy = RetryPolicy {
        max_retries: 3,
        initial_delay: Duration::from_millis(50),
        max_delay: Duration::from_millis(200),
        per_attempt_timeout: None,
    };
    let mut client = PredictClient::connect_tls_with_retry(endpoint, &client_tls, policy)
        .await
        .expect("connect_tls");

    let req = make_request_i8("r-mtls", "echo-model", vec![1, 1], &[1]);
    let resp = client.predict(req).await.expect("predict");

    assert_eq!(resp.request_id, "r-mtls");
    assert_eq!(resp.status, PredictStatus::Ok as i32);

    server.abort();
}

#[tokio::test]
async fn client_with_wrong_ca_fails_connect() {
    let certs = gen_ca_and_identities();
    let server_tls = ServerTlsOptions {
        cert_pem: certs.server_cert_pem.clone(),
        key_pem: certs.server_key_pem.clone(),
        client_ca_pem: Some(certs.ca_pem.clone()),
    };
    let (endpoint, server) = spawn_mtls_server(server_tls).await;

    let client_tls = ClientTlsOptions {
        server_ca_pem: certs.other_ca_pem.clone(),
        client_identity_pem: Some((
            certs.client_cert_pem.clone(),
            certs.client_key_pem.clone(),
        )),
        domain: "localhost".to_string(),
    };
    let res = PredictClient::connect_tls_with_retry(
        endpoint,
        &client_tls,
        RetryPolicy::none(),
    )
    .await;
    let err = res.err().expect("must fail");
    assert!(
        matches!(err, PredictError::Transport(_)),
        "expected Transport error, got {}",
        err
    );

    server.abort();
}

#[tokio::test]
async fn client_without_identity_rejected_under_mtls() {
    let certs = gen_ca_and_identities();
    let server_tls = ServerTlsOptions {
        cert_pem: certs.server_cert_pem.clone(),
        key_pem: certs.server_key_pem.clone(),
        client_ca_pem: Some(certs.ca_pem.clone()),
    };
    let (endpoint, server) = spawn_mtls_server(server_tls).await;

    let client_tls = ClientTlsOptions {
        server_ca_pem: certs.ca_pem.clone(),
        client_identity_pem: None,
        domain: "localhost".to_string(),
    };
    let res = PredictClient::connect_tls_with_retry(
        endpoint.clone(),
        &client_tls,
        RetryPolicy::none(),
    )
    .await;

    match res {
        Ok(mut client) => {
            let req = make_request_i8("noid", "echo-model", vec![1, 1], &[1]);
            let err = client
                .predict(req)
                .await
                .expect_err("predict must fail under mTLS without client cert");
            match err {
                PredictError::Transport(_) => {  }
                PredictError::Grpc(status) => {
                    let code = status.code();
                    assert!(
                        matches!(
                            code,
                            tonic::Code::Unavailable
                                | tonic::Code::Unknown
                                | tonic::Code::Internal
                                | tonic::Code::Cancelled
                        ),
                        "expected handshake-class gRPC code, got {:?}: {}",
                        code,
                        status.message()
                    );
                }
                other => panic!("unexpected error type: {}", other),
            }
        }
        Err(PredictError::Transport(_)) => {  }
        Err(other) => panic!("unexpected error {}", other),
    }

    server.abort();
}

#[tokio::test]
async fn client_with_domain_mismatch_fails_connect() {

    let certs = gen_ca_and_identities();
    let server_tls = ServerTlsOptions {
        cert_pem: certs.server_cert_pem.clone(),
        key_pem: certs.server_key_pem.clone(),
        client_ca_pem: Some(certs.ca_pem.clone()),
    };
    let (endpoint, server) = spawn_mtls_server(server_tls).await;

    let client_tls = ClientTlsOptions {
        server_ca_pem: certs.ca_pem.clone(),
        client_identity_pem: Some((
            certs.client_cert_pem.clone(),
            certs.client_key_pem.clone(),
        )),
        domain: "not-localhost.example".to_string(),
    };
    let res = PredictClient::connect_tls_with_retry(
        endpoint,
        &client_tls,
        RetryPolicy::none(),
    )
    .await;
    let err = res.err().expect("must fail on hostname mismatch");
    assert!(
        matches!(err, PredictError::Transport(_)),
        "expected Transport error from hostname verification, got {}",
        err
    );

    server.abort();
}
