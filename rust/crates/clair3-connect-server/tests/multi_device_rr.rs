use clair3_connect_infer::model_registry::Registry;
use clair3_connect_server::metrics_server::MetricsSnapshot;
use clair3_connect_server::service::PredictService;
use clair3_connect_transport::codec::decode_f32;
use clair3_connect_transport::predict_client::{make_request_i32, PredictClient};
use clair3_connect_transport::proto::clair3_connect_inference_server::Clair3ConnectInferenceServer;
use std::path::PathBuf;
use std::time::Duration;
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

#[tokio::test(flavor = "multi_thread", worker_threads = 4)]
async fn round_robin_distributes_load_across_two_workers() {
    let yaml = workspace_root().join("models/model_registry.example.yaml");
    if !yaml.exists() {
        eprintln!("Skipping — registry yaml missing");
        return;
    }
    let reg = Registry::load(&yaml).expect("registry");
    let yaml_dir = yaml.parent().unwrap().to_path_buf();

    let svc = match PredictService::load_with_options(
        &reg,
        &yaml_dir,
        &[tch::Device::Cpu, tch::Device::Cpu],
        clair3_connect_server::backpressure::BackpressureConfig {
            max_inflight: 8,
            retry_after: Duration::from_millis(100),
        },
        clair3_connect_infer::batcher::BatchConfig::default(),
    ) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Skipping — model load failed: {e}");
            return;
        }
    };

    let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    let incoming = tokio_stream::wrappers::TcpListenerStream::new(listener);
    let svc_for_server = svc.clone();
    let server = tokio::spawn(async move {
        Server::builder()
            .add_service(Clair3ConnectInferenceServer::new(svc_for_server))
            .serve_with_incoming(incoming)
            .await
            .ok();
    });
    tokio::task::yield_now().await;

    let endpoint = format!("http://{}", addr);
    let n: usize = 16;
    let mut handles = Vec::with_capacity(n);
    for i in 0..n {
        let endpoint = endpoint.clone();
        handles.push(tokio::spawn(async move {
            let mut client = PredictClient::connect(endpoint).await.expect("connect");
            let zeros = vec![0_i32; 33 * 18];
            let req = make_request_i32(
                &format!("rr{i}"),
                "pileup_ont_hac_v500",
                vec![1, 33, 18],
                &zeros,
            );
            let resp = client.predict(req).await.expect("predict");
            let y = decode_f32(&resp.y_payload);
            assert_eq!(y.len(), 24);
        }));
    }
    for h in handles {
        h.await.unwrap();
    }

    let snap: MetricsSnapshot = svc.metrics_snapshot();
    let pileup_rows: Vec<_> = snap
        .per_model
        .iter()
        .filter(|m| m.model_id == "pileup_ont_hac_v500")
        .collect();
    assert_eq!(
        pileup_rows.len(),
        2,
        "expected 2 entries (one per device); got {}",
        pileup_rows.len()
    );
    let total_items: u64 = pileup_rows.iter().map(|m| m.items).sum();
    assert_eq!(total_items, n as u64, "total items must equal request count");
    for m in &pileup_rows {
        assert!(
            m.items > 0,
            "worker {} got no items — round-robin broken",
            m.device
        );
        eprintln!(
            "device={} batches={} items={}",
            m.device, m.batches, m.items
        );
    }

    server.abort();
}
