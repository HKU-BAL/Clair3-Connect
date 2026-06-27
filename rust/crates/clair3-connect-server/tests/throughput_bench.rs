use clair3_connect_infer::batcher::BatchConfig;
use clair3_connect_infer::model_registry::Registry;
use clair3_connect_server::backpressure::BackpressureConfig;
use clair3_connect_server::service::PredictService;
use clair3_connect_transport::codec::decode_f32;
use clair3_connect_transport::predict_client::{make_request_i32, PredictClient};
use clair3_connect_transport::proto::clair3_connect_inference_server::Clair3ConnectInferenceServer;
use std::path::PathBuf;
use std::time::{Duration, Instant};
use tokio::net::TcpListener;
use tonic::transport::Server;

const N_CLIENTS: usize = 8;
const REQUESTS_PER_CLIENT: usize = 32;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

async fn run_one(svc: PredictService) -> (usize, Duration) {
    let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    let incoming = tokio_stream::wrappers::TcpListenerStream::new(listener);
    let server = tokio::spawn(async move {
        Server::builder()
            .add_service(Clair3ConnectInferenceServer::new(svc))
            .serve_with_incoming(incoming)
            .await
            .ok();
    });

    tokio::task::yield_now().await;

    let endpoint = format!("http://{}", addr);
    let total = N_CLIENTS * REQUESTS_PER_CLIENT;
    let start = Instant::now();
    let mut handles = Vec::with_capacity(N_CLIENTS);
    for client_idx in 0..N_CLIENTS {
        let endpoint = endpoint.clone();
        handles.push(tokio::spawn(async move {
            let mut client = PredictClient::connect(endpoint).await.expect("connect");
            for r in 0..REQUESTS_PER_CLIENT {

                let zeros = vec![0_i32; 33 * 18];
                let req = make_request_i32(
                    &format!("c{client_idx}_r{r}"),
                    "pileup_ont_hac_v500",
                    vec![1, 33, 18],
                    &zeros,
                );
                let resp = client.predict(req).await.expect("predict");
                let y = decode_f32(&resp.y_payload);
                assert_eq!(y.len(), 24, "expected 1 row x 24 cols");
            }
        }));
    }
    for h in handles {
        h.await.unwrap();
    }
    let elapsed = start.elapsed();
    server.abort();
    (total, elapsed)
}

#[tokio::test(flavor = "multi_thread", worker_threads = 8)]
#[ignore = "benchmark — runs ~10s, invoke explicitly"]
async fn batching_beats_no_batching_by_at_least_1_5x() {
    let yaml = workspace_root().join("models/model_registry.example.yaml");
    if !yaml.exists() {
        eprintln!("Skipping — registry yaml missing");
        return;
    }
    let reg = Registry::load(&yaml).expect("registry");
    let yaml_dir = yaml.parent().unwrap().to_path_buf();

    let baseline_bp = BackpressureConfig {
        max_inflight: N_CLIENTS * 2,
        retry_after: Duration::from_millis(100),
    };
    let baseline_batch = BatchConfig {
        max_batch_size: 1,
        flush_timeout: Duration::from_millis(1),
        queue_depth: 1024,
    };
    let baseline_svc = match PredictService::load_with_options(
        &reg,
        &yaml_dir,
        &[tch::Device::Cpu],
        baseline_bp.clone(),
        baseline_batch,
    ) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Skipping — model load failed: {e}");
            return;
        }
    };
    let (baseline_total, baseline_elapsed) = run_one(baseline_svc).await;
    let baseline_rps = baseline_total as f64 / baseline_elapsed.as_secs_f64();
    eprintln!(
        "BASELINE (no batching): {} requests in {:?} = {:.1} rps",
        baseline_total, baseline_elapsed, baseline_rps
    );

    let batched_bp = baseline_bp.clone();
    let batched_batch = BatchConfig {
        max_batch_size: 64,
        flush_timeout: Duration::from_millis(4),
        queue_depth: 1024,
    };
    let batched_svc = PredictService::load_with_options(
        &reg,
        &yaml_dir,
        &[tch::Device::Cpu],
        batched_bp,
        batched_batch,
    )
    .expect("load batched");
    let (batched_total, batched_elapsed) = run_one(batched_svc).await;
    let batched_rps = batched_total as f64 / batched_elapsed.as_secs_f64();
    eprintln!(
        "BATCHED (max_batch=64, flush=4ms): {} requests in {:?} = {:.1} rps",
        batched_total, batched_elapsed, batched_rps
    );

    let speedup = batched_rps / baseline_rps;
    eprintln!("speedup = {:.2}×", speedup);

    assert!(
        speedup >= 1.5,
        "expected ≥1.5× speedup from batching; got {:.2}× (baseline {:.1} rps vs batched {:.1} rps)",
        speedup,
        baseline_rps,
        batched_rps
    );
}
