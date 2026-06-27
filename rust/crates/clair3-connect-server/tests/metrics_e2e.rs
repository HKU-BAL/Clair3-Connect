use clair3_connect_infer::model_registry::Registry;
use clair3_connect_server::metrics_server;
use clair3_connect_server::service::PredictService;
use clair3_connect_transport::codec::decode_f32;
use clair3_connect_transport::predict_client::{make_request_i32, PredictClient};
use clair3_connect_transport::proto::clair3_connect_inference_server::Clair3ConnectInferenceServer;
use std::path::PathBuf;
use tokio::io::{AsyncReadExt, AsyncWriteExt};
use tokio::net::{TcpListener, TcpStream};
use tonic::transport::Server;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

async fn http_get(addr: std::net::SocketAddr, path: &str) -> (String, u16) {
    let mut s = TcpStream::connect(addr).await.unwrap();
    let req = format!("GET {path} HTTP/1.1\r\nHost: localhost\r\nConnection: close\r\n\r\n");
    s.write_all(req.as_bytes()).await.unwrap();
    let mut buf = Vec::new();
    s.read_to_end(&mut buf).await.unwrap();
    let text = String::from_utf8_lossy(&buf).to_string();

    let code: u16 = text
        .lines()
        .next()
        .and_then(|l| l.split_whitespace().nth(1))
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);

    let body = text
        .split_once("\r\n\r\n")
        .map(|(_, b)| b.to_string())
        .unwrap_or_default();
    (body, code)
}

#[tokio::test]
async fn metrics_endpoint_reports_predict_activity() {
    let yaml = workspace_root().join("models/model_registry.example.yaml");
    let x_path =
        workspace_root().join("test_data/ont_chr20_demo/artifacts/pileup_post_scale_X_int32.npy");
    if !yaml.exists() || !x_path.exists() {
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

    let metrics_listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let metrics_addr = metrics_listener.local_addr().unwrap();
    drop(metrics_listener);
    let _metrics_handle = metrics_server::spawn(metrics_addr, svc.clone())
        .await
        .expect("metrics_server::spawn");

    let (body, code) = http_get(metrics_addr, "/metrics").await;
    assert_eq!(code, 200);
    assert!(body.contains("# TYPE clair3_connect_batches_total counter"));
    assert!(
        body.contains("clair3_connect_batches_total{model_id=\"pileup_ont_hac_v500\""),
        "expected pileup model series; got:\n{body}"
    );

    for line in body.lines() {
        if line.starts_with("clair3_connect_batches_total{") {
            assert!(
                line.ends_with(" 0"),
                "expected batches=0 pre-request; got: {line}"
            );
        }
    }

    let (body, code) = http_get(metrics_addr, "/healthz").await;
    assert_eq!(code, 200);
    assert_eq!(body, "ok\n");

    let (_body, code) = http_get(metrics_addr, "/does/not/exist").await;
    assert_eq!(code, 404);

    let grpc_listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let grpc_addr = grpc_listener.local_addr().unwrap();
    let incoming = tokio_stream::wrappers::TcpListenerStream::new(grpc_listener);
    let server = tokio::spawn(async move {
        Server::builder()
            .add_service(Clair3ConnectInferenceServer::new(svc))
            .serve_with_incoming(incoming)
            .await
            .ok();
    });

    let endpoint = format!("http://{}", grpc_addr);
    let mut client = PredictClient::connect(endpoint).await.expect("connect");
    let x_bytes = std::fs::read(&x_path).unwrap();

    let major = x_bytes[6];
    let header_len = if major == 1 {
        u16::from_le_bytes([x_bytes[8], x_bytes[9]]) as usize + 10
    } else {
        u32::from_le_bytes([x_bytes[8], x_bytes[9], x_bytes[10], x_bytes[11]]) as usize + 12
    };
    let x_int32: Vec<i32> = x_bytes[header_len..]
        .chunks_exact(4)
        .map(|c| i32::from_le_bytes([c[0], c[1], c[2], c[3]]))
        .take(33 * 18)
        .collect();
    let req = make_request_i32("m1", "pileup_ont_hac_v500", vec![1, 33, 18], &x_int32);
    let resp = client.predict(req).await.expect("predict");
    let y = decode_f32(&resp.y_payload);
    assert_eq!(y.len(), 24, "1 row x 24 cols");

    let (body, _code) = http_get(metrics_addr, "/metrics").await;
    let line = body
        .lines()
        .find(|l| l.starts_with("clair3_connect_items_total{model_id=\"pileup_ont_hac_v500\""))
        .expect("items_total series present");
    let value: u64 = line
        .rsplit_once(' ')
        .map(|(_, v)| v.parse().unwrap_or(0))
        .unwrap_or(0);
    assert!(
        value >= 1,
        "expected items_total >= 1 after Predict; got line: {line}"
    );

    server.abort();
}
