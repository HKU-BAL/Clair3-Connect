use clair3_connect_infer::model_registry::Registry;
use clair3_connect_server::service::PredictService;
use clair3_connect_transport::codec::decode_f32;
use clair3_connect_transport::predict_client::{make_request_i32, PredictClient};
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
    assert_eq!(&bytes[0..6], b"\x93NUMPY");
    let major = bytes[6];
    let header_len = if major == 1 {
        u16::from_le_bytes([bytes[8], bytes[9]]) as usize + 10
    } else {
        u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize + 12
    };
    let data = &bytes[header_len..];
    data.chunks_exact(4)
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
    let data = &bytes[header_len..];
    data.chunks_exact(4)
        .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .collect()
}

#[tokio::test]
async fn predict_roundtrip_pileup() {
    let yaml = workspace_root().join("models/model_registry.example.yaml");
    let x_path = workspace_root().join("test_data/ont_chr20_demo/artifacts/pileup_post_scale_X_int32.npy");
    let y_path = workspace_root().join("test_data/ont_chr20_demo/artifacts/Y_pileup_cpu.npy");
    if !yaml.exists() || !x_path.exists() || !y_path.exists() {
        eprintln!("Skipping — fixtures missing");
        return;
    }

    let reg = Registry::load(&yaml).expect("registry");
    let yaml_dir = yaml.parent().unwrap().to_path_buf();
    let svc = match PredictService::load(&reg, &yaml_dir, tch::Device::Cpu) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Skipping — model load failed: {e}");
            return;
        }
    };

    let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    let incoming = tokio_stream::wrappers::TcpListenerStream::new(listener);

    let server_handle = tokio::spawn(async move {
        Server::builder()
            .add_service(Clair3ConnectInferenceServer::new(svc))
            .serve_with_incoming(incoming)
            .await
            .ok();
    });

    let endpoint = format!("http://{}", addr);
    let mut client = PredictClient::connect(endpoint).await.expect("connect");

    let x_int32 = read_npy_i32(&x_path);
    let n: i64 = 988;
    let d1: i64 = 33;
    let d2: i64 = 18;
    let elements_per_row: i64 = d1 * d2;
    assert_eq!(x_int32.len() as i64, n * elements_per_row, "X size unexpected");

    let batch_size: i64 = 32;
    let mut y_all: Vec<f32> = Vec::with_capacity((n * 24) as usize);
    let mut row = 0_i64;
    while row < n {
        let end = (row + batch_size).min(n);
        let batch = end - row;
        let start_byte = (row * elements_per_row) as usize;
        let end_byte = (end * elements_per_row) as usize;
        let slice = &x_int32[start_byte..end_byte];
        let req = make_request_i32(
            &format!("batch{}", row),
            "pileup_ont_hac_v500",
            vec![batch, d1, d2],
            slice,
        );
        let resp = client.predict(req).await.expect("predict");
        let y_batch = decode_f32(&resp.y_payload);
        assert_eq!(y_batch.len() as i64, batch * 24, "Y batch shape unexpected");
        y_all.extend_from_slice(&y_batch);
        row = end;
    }

    let y_py = read_npy_f32(&y_path);
    assert_eq!(y_all.len(), y_py.len(), "Y total len mismatch");

    let mut max_abs = 0.0_f32;
    for (a, b) in y_all.iter().zip(y_py.iter()) {
        let d = (a - b).abs();
        if d > max_abs {
            max_abs = d;
        }
    }
    eprintln!("e2e max_abs_diff = {max_abs:e}");
    assert!(max_abs < 1e-5, "predict e2e parity exceeded: {max_abs:e}");

    server_handle.abort();
}
