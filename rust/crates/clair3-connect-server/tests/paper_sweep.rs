use clair3_connect_infer::batcher::BatchConfig;
use clair3_connect_infer::model_registry::Registry;
use clair3_connect_server::backpressure::BackpressureConfig;
use clair3_connect_server::service::PredictService;
use clair3_connect_transport::codec::decode_f32;
use clair3_connect_transport::predict_client::{make_request_i32, make_request_i8, PredictClient};
use clair3_connect_transport::proto::clair3_connect_inference_server::Clair3ConnectInferenceServer;
use std::env;
use std::path::PathBuf;
use std::time::{Duration, Instant};
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

fn parse_devices(s: &str) -> Vec<tch::Device> {
    s.split(',')
        .map(|t| t.trim())
        .filter(|t| !t.is_empty())
        .map(|t| {
            if t == "cpu" {
                tch::Device::Cpu
            } else if let Some(n) = t.strip_prefix("cuda:") {
                tch::Device::Cuda(n.parse().expect("cuda index"))
            } else {
                panic!("bad device token: {t}")
            }
        })
        .collect()
}

fn envv(k: &str, d: &str) -> String {
    env::var(k).unwrap_or_else(|_| d.to_string())
}

async fn run_load(
    svc: PredictService,
    n_clients: usize,
    reqs: usize,
    model_id: String,
) -> (usize, Duration, Vec<u64>) {
    let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    let incoming = tokio_stream::wrappers::TcpListenerStream::new(listener);
    let server = tokio::spawn(async move {
        Server::builder()
            .tcp_nodelay(true)
            .add_service(Clair3ConnectInferenceServer::new(svc))
            .serve_with_incoming(incoming)
            .await
            .ok();
    });
    tokio::task::yield_now().await;

    let endpoint = format!("http://{}", addr);
    let total = n_clients * reqs;
    let start = Instant::now();
    let mut handles = Vec::with_capacity(n_clients);
    for client_idx in 0..n_clients {
        let endpoint = endpoint.clone();
        let model_id = model_id.clone();
        handles.push(tokio::spawn(async move {
            let mut client = PredictClient::connect(endpoint).await.expect("connect");
            let mut lats = Vec::with_capacity(reqs);

            let is_fa = model_id.contains("full_alignment");
            for r in 0..reqs {
                let rid = format!("c{client_idx}_r{r}");
                let (req, expect_len) = if is_fa {
                    let zeros = vec![0_i8; 89 * 33 * 8];
                    (make_request_i8(&rid, &model_id, vec![1, 89, 33, 8], &zeros), 90)
                } else {
                    let zeros = vec![0_i32; 33 * 18];
                    (make_request_i32(&rid, &model_id, vec![1, 33, 18], &zeros), 24)
                };
                let t0 = Instant::now();
                let resp = client.predict(req).await.expect("predict");
                let dt = t0.elapsed().as_micros() as u64;
                let y = decode_f32(&resp.y_payload);
                assert_eq!(y.len(), expect_len);
                lats.push(dt);
            }
            lats
        }));
    }
    let mut all_lats = Vec::with_capacity(total);
    for h in handles {
        all_lats.extend(h.await.unwrap());
    }
    let elapsed = start.elapsed();
    server.abort();
    (total, elapsed, all_lats)
}

fn pct(sorted: &[u64], p: f64) -> u64 {
    if sorted.is_empty() {
        return 0;
    }
    let idx = ((p / 100.0) * (sorted.len() as f64 - 1.0)).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

#[tokio::test(flavor = "multi_thread", worker_threads = 16)]
#[ignore = "paper benchmark — invoke explicitly with env knobs"]
async fn paper_throughput_latency_sweep() {
    let device_str = envv("CLAIR3_CONNECT_BENCH_DEVICE", "cpu");
    let devices = parse_devices(&device_str);
    let out_path = envv("CLAIR3_CONNECT_BENCH_OUT", "/tmp/paper_sweep.json");
    let n_clients: usize = envv("CLAIR3_CONNECT_BENCH_CLIENTS", "16").parse().unwrap();
    let reqs: usize = envv("CLAIR3_CONNECT_BENCH_REQS", "64").parse().unwrap();
    let model_id = envv("CLAIR3_CONNECT_BENCH_MODEL", "pileup_ont_hac_v500");
    let batches: Vec<usize> = envv("CLAIR3_CONNECT_BENCH_BATCHES", "1,2,4,8,16,32,64")
        .split(',')
        .map(|x| x.trim().parse().unwrap())
        .collect();
    let registry = envv(
        "CLAIR3_CONNECT_BENCH_REGISTRY",
        workspace_root()
            .join("models/model_registry.example.yaml")
            .to_str()
            .unwrap(),
    );

    let yaml = PathBuf::from(&registry);
    let reg = Registry::load(&yaml).expect("registry");
    let yaml_dir = yaml.parent().unwrap().to_path_buf();

    eprintln!(
        "paper_sweep: device={device_str} clients={n_clients} reqs={reqs} model={model_id} batches={batches:?}"
    );

    {
        let bp = BackpressureConfig {
            max_inflight: n_clients * 2,
            retry_after: Duration::from_millis(100),
        };
        let warm = BatchConfig {
            max_batch_size: 16,
            flush_timeout: Duration::from_millis(5),
            queue_depth: 1024,
        };
        if let Ok(svc) =
            PredictService::load_with_options(&reg, &yaml_dir, &devices, bp, warm)
        {
            let _ = run_load(svc, n_clients, 8.min(reqs).max(1), model_id.clone()).await;
            eprintln!("warmup done");
        } else {
            eprintln!("Skipping — model load failed (device {device_str} unavailable?)");
            return;
        }
    }

    let mut records: Vec<String> = Vec::new();
    for &mb in &batches {
        let bp = BackpressureConfig {
            max_inflight: n_clients * 2,
            retry_after: Duration::from_millis(100),
        };
        let batch = BatchConfig {
            max_batch_size: mb,

            flush_timeout: if mb == 1 {
                Duration::from_nanos(1)
            } else {
                Duration::from_millis(5)
            },
            queue_depth: 1024,
        };
        let svc =
            match PredictService::load_with_options(&reg, &yaml_dir, &devices, bp, batch) {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("Skipping mb={mb} — load failed: {e}");
                    continue;
                }
            };
        let (total, elapsed, mut lats) = run_load(svc, n_clients, reqs, model_id.clone()).await;
        let rps = total as f64 / elapsed.as_secs_f64();
        lats.sort_unstable();
        let mean = lats.iter().sum::<u64>() as f64 / lats.len().max(1) as f64;
        eprintln!(
            "mb={mb:>3}  rps={rps:>9.1}  p50={}us p95={}us p99={}us",
            pct(&lats, 50.0),
            pct(&lats, 95.0),
            pct(&lats, 99.0)
        );
        records.push(format!(
            "{{\"max_batch\":{mb},\"clients\":{n_clients},\"reqs_per_client\":{reqs},\
\"total_requests\":{total},\"elapsed_s\":{:.6},\"rps\":{:.3},\
\"lat_us_mean\":{:.1},\"lat_us_p50\":{},\"lat_us_p95\":{},\"lat_us_p99\":{},\
\"device\":\"{}\",\"model\":\"{}\"}}",
            elapsed.as_secs_f64(),
            rps,
            mean,
            pct(&lats, 50.0),
            pct(&lats, 95.0),
            pct(&lats, 99.0),
            device_str,
            model_id,
        ));
    }

    let json = format!("[\n  {}\n]\n", records.join(",\n  "));
    std::fs::write(&out_path, json).expect("write json");
    eprintln!("wrote {out_path}");
}
