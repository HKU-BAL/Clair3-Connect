use clair3_connect_infer::model_registry::Registry;
use clair3_connect_server::service::PredictService;
use clair3_connect_transport::aead::AeadKey;
use clair3_connect_transport::codec::decode_f32;
use clair3_connect_transport::predict_client::{
    make_request_i32, PredictClient, PredictError,
};
use clair3_connect_transport::proto::clair3_connect_inference_server::Clair3ConnectInferenceServer;
use std::path::PathBuf;
use tokio::net::TcpListener;
use tonic::transport::Server;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

fn read_npy_i32(path: &std::path::Path) -> Vec<i32> {
    let bytes = std::fs::read(path).expect("read npy");
    let major = bytes[6];
    let header_len = if major == 1 {
        u16::from_le_bytes([bytes[8], bytes[9]]) as usize + 10
    } else {
        u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize + 12
    };
    bytes[header_len..]
        .chunks_exact(4)
        .map(|c| i32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect()
}

fn read_npy_f32(path: &std::path::Path) -> Vec<f32> {
    let bytes = std::fs::read(path).expect("read npy");
    let major = bytes[6];
    let header_len = if major == 1 {
        u16::from_le_bytes([bytes[8], bytes[9]]) as usize + 10
    } else {
        u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize + 12
    };
    bytes[header_len..]
        .chunks_exact(4)
        .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect()
}

fn load_service_or_skip(server_key: Option<AeadKey>) -> Option<PredictService> {
    let yaml = workspace_root().join("models/model_registry.example.yaml");
    if !yaml.exists() {
        eprintln!("Skipping — registry yaml missing");
        return None;
    }
    let reg = Registry::load(&yaml).ok()?;
    let yaml_dir = yaml.parent()?.to_path_buf();
    let mut svc = PredictService::load(&reg, &yaml_dir, tch::Device::Cpu).ok()?;
    if let Some(k) = server_key {
        svc.set_aead_key(k);
    }
    Some(svc)
}

async fn spawn_grpc_server(svc: PredictService) -> (String, tokio::task::JoinHandle<()>) {
    let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    let incoming = tokio_stream::wrappers::TcpListenerStream::new(listener);
    let handle = tokio::spawn(async move {
        Server::builder()
            .add_service(Clair3ConnectInferenceServer::new(svc))
            .serve_with_incoming(incoming)
            .await
            .ok();
    });
    tokio::task::yield_now().await;
    (format!("http://{}", addr), handle)
}

#[tokio::test]
async fn aead_predict_roundtrip() {

    let key_bytes = [0xa5u8; 32];
    let server_key = AeadKey::from_bytes(&key_bytes).unwrap();
    let client_key = AeadKey::from_bytes(&key_bytes).unwrap();

    let svc = match load_service_or_skip(Some(server_key)) {
        Some(s) => s,
        None => return,
    };
    let (endpoint, server) = spawn_grpc_server(svc).await;
    let mut client = PredictClient::connect(endpoint).await.expect("connect");
    client.set_aead_key(client_key);
    assert!(client.aead_enabled());

    let x_path = workspace_root().join("test_data/ont_chr20_demo/artifacts/pileup_post_scale_X_int32.npy");
    let y_path = workspace_root().join("test_data/ont_chr20_demo/artifacts/Y_pileup_cpu.npy");
    if !x_path.exists() || !y_path.exists() {
        eprintln!("Skipping — fixture .npy missing");
        server.abort();
        return;
    }
    let x_int32 = read_npy_i32(&x_path);
    let row: Vec<i32> = x_int32[..33 * 18].to_vec();
    let req = make_request_i32("aead-rt", "pileup_ont_hac_v500", vec![1, 33, 18], &row);
    let resp = client.predict(req).await.expect("predict");
    let y = decode_f32(&resp.y_payload);
    let y_py = read_npy_f32(&y_path);
    let mut max_diff = 0.0_f32;
    for i in 0..24 {
        max_diff = max_diff.max((y[i] - y_py[i]).abs());
    }
    eprintln!("aead_e2e max_abs_diff over first row = {max_diff:e}");
    assert!(max_diff < 1e-5, "AEAD round-trip should preserve tensor bytes");

    server.abort();
}

#[tokio::test]
async fn aead_required_server_rejects_plaintext() {
    let key_bytes = [0xb6u8; 32];
    let server_key = AeadKey::from_bytes(&key_bytes).unwrap();
    let svc = match load_service_or_skip(Some(server_key)) {
        Some(s) => s,
        None => return,
    };
    let (endpoint, server) = spawn_grpc_server(svc).await;

    let mut client = PredictClient::connect(endpoint).await.expect("connect");
    let req = make_request_i32(
        "no-key",
        "pileup_ont_hac_v500",
        vec![1, 33, 18],
        &vec![0_i32; 33 * 18],
    );
    let err = client.predict(req).await.expect_err("must reject");
    match err {
        PredictError::NonOk { code, .. } => {
            assert_eq!(code, "aead_required", "wrong error_code");
        }
        other => panic!("expected NonOk(aead_required), got {other}"),
    }
    server.abort();
}

#[tokio::test]
async fn aead_not_supported_server_rejects_ciphertext() {

    let svc = match load_service_or_skip(None) {
        Some(s) => s,
        None => return,
    };
    let (endpoint, server) = spawn_grpc_server(svc).await;
    let mut client = PredictClient::connect(endpoint).await.expect("connect");
    let client_key = AeadKey::from_bytes(&[0xc7u8; 32]).unwrap();
    client.set_aead_key(client_key);
    let req = make_request_i32(
        "ct",
        "pileup_ont_hac_v500",
        vec![1, 33, 18],
        &vec![0_i32; 33 * 18],
    );
    let err = client.predict(req).await.expect_err("must reject");
    match err {
        PredictError::NonOk { code, .. } => {
            assert_eq!(code, "aead_not_supported", "wrong error_code");
        }
        other => panic!("expected NonOk(aead_not_supported), got {other}"),
    }
    server.abort();
}

#[tokio::test]
async fn aead_wrong_key_rejected() {
    let server_key = AeadKey::from_bytes(&[0xd8u8; 32]).unwrap();
    let client_key = AeadKey::from_bytes(&[0xe9u8; 32]).unwrap();
    let svc = match load_service_or_skip(Some(server_key)) {
        Some(s) => s,
        None => return,
    };
    let (endpoint, server) = spawn_grpc_server(svc).await;
    let mut client = PredictClient::connect(endpoint).await.expect("connect");
    client.set_aead_key(client_key);
    let req = make_request_i32(
        "wrong-key",
        "pileup_ont_hac_v500",
        vec![1, 33, 18],
        &vec![0_i32; 33 * 18],
    );
    let err = client.predict(req).await.expect_err("must reject");
    match err {
        PredictError::NonOk { code, .. } => {
            assert_eq!(code, "aead_open_failed", "wrong error_code");
        }
        other => panic!("expected NonOk(aead_open_failed), got {other}"),
    }
    server.abort();
}
