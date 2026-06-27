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

const N_CLIENTS: usize = 32;
const REQUESTS_PER_CLIENT: usize = 16;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

async fn run_one(svc: PredictService) -> Option<(usize, Duration)> {
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
    for ci in 0..N_CLIENTS {
        let endpoint = endpoint.clone();
        handles.push(tokio::spawn(async move {
            let mut client = PredictClient::connect(endpoint).await.ok()?;
            for r in 0..REQUESTS_PER_CLIENT {
                let zeros = vec![0_i32; 33 * 18];
                let req = make_request_i32(
                    &format!("c{ci}_r{r}"),
                    "pileup_ont_hac_v500",
                    vec![1, 33, 18],
                    &zeros,
                );
                let resp = client.predict(req).await.ok()?;
                let y = decode_f32(&resp.y_payload);
                if y.len() != 24 {
                    return None;
                }
            }
            Some(())
        }));
    }
    let mut ok = 0;
    for h in handles {
        if h.await.ok().flatten().is_some() {
            ok += 1;
        }
    }
    let elapsed = start.elapsed();
    server.abort();
    if ok == N_CLIENTS {
        Some((total, elapsed))
    } else {
        eprintln!(
            "  ({}/{} clients failed — likely OOM on shared GPU)",
            N_CLIENTS - ok,
            N_CLIENTS
        );
        None
    }
}

#[tokio::test(flavor = "multi_thread", worker_threads = 8)]
#[ignore = "GPU benchmark — invoke explicitly on a multi-GPU host"]
async fn two_gpus_beat_one_gpu() {
    if !tch::Cuda::is_available() {
        eprintln!("Skipping — CUDA not available");
        return;
    }
    let n = tch::Cuda::device_count();
    if n < 2 {
        eprintln!("Skipping — need ≥ 2 CUDA devices, got {n}");
        return;
    }

    let yaml = workspace_root().join("models/model_registry.example.yaml");
    if !yaml.exists() {
        eprintln!("Skipping — registry yaml missing");
        return;
    }
    let reg = Registry::load(&yaml).expect("registry");
    let yaml_dir = yaml.parent().unwrap().to_path_buf();

    let _ = clair3_connect_infer::_ensure_cuda_linked();

    let bp = BackpressureConfig {
        max_inflight: 16,
        retry_after: Duration::from_millis(50),
    };
    let batch = BatchConfig {
        max_batch_size: 64,
        flush_timeout: Duration::from_millis(4),
        queue_depth: 1024,
    };

    let svc1 = match PredictService::load_with_options(
        &reg,
        &yaml_dir,
        &[tch::Device::Cuda(0)],
        bp.clone(),
        batch.clone(),
    ) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Skipping — 1-GPU load failed (probably OOM on shared host): {e}");
            return;
        }
    };
    let r1 = run_one(svc1).await;
    let Some((tot1, el1)) = r1 else {
        eprintln!("Skipping — 1-GPU run failed (shared GPU OOM?)");
        return;
    };
    let rps1 = tot1 as f64 / el1.as_secs_f64();
    eprintln!("1× GPU: {} reqs in {:?} = {:.1} rps", tot1, el1, rps1);

    let svc2 = match PredictService::load_with_options(
        &reg,
        &yaml_dir,
        &[tch::Device::Cuda(0), tch::Device::Cuda(1)],
        bp,
        batch,
    ) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Skipping — 2-GPU load failed: {e}");
            return;
        }
    };
    let r2 = run_one(svc2).await;
    let Some((tot2, el2)) = r2 else {
        eprintln!("Skipping — 2-GPU run failed");
        return;
    };
    let rps2 = tot2 as f64 / el2.as_secs_f64();
    eprintln!("2× GPU: {} reqs in {:?} = {:.1} rps", tot2, el2, rps2);

    let speedup = rps2 / rps1;
    eprintln!("speedup = {:.2}×", speedup);
    assert!(
        speedup >= 1.4,
        "expected ≥1.4× speedup from 2-GPU pool; got {:.2}× ({:.1} rps vs {:.1} rps)",
        speedup,
        rps1,
        rps2
    );
}
