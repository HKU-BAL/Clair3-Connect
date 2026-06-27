use clair3_connect_transport::aead::AeadKey;
use clair3_connect_transport::predict_client::{make_request_i8, PredictClient, RetryPolicy};
use clair3_connect_transport::tls::ClientTlsOptions;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

fn load_aead() -> Option<AeadKey> {
    let p = std::env::var("FA_AEAD_KEY").ok()?;
    let bytes = std::fs::read(&p).expect("read FA_AEAD_KEY");
    Some(AeadKey::from_bytes(&bytes).expect("FA_AEAD_KEY must be 32 bytes"))
}

fn tls_opts() -> Option<ClientTlsOptions> {
    let ca = std::env::var("FA_TLS_CA").ok()?;
    let domain = std::env::var("FA_TLS_DOMAIN").unwrap_or_else(|_| "localhost".into());
    let cert = std::env::var("FA_TLS_CERT").ok();
    let key = std::env::var("FA_TLS_KEY").ok();
    Some(
        ClientTlsOptions::from_paths(
            Path::new(&ca),
            cert.as_deref().map(Path::new),
            key.as_deref().map(Path::new),
            &domain,
        )
        .expect("build ClientTlsOptions"),
    )
}

async fn make_client(endpoint: &str) -> PredictClient {
    let mut c = match tls_opts() {
        Some(tls) => {
            PredictClient::connect_tls_with_retry(endpoint.to_string(), &tls, RetryPolicy::none())
                .await
                .expect("connect_tls")
        }
        None => PredictClient::connect(endpoint.to_string())
            .await
            .expect("connect"),
    };
    if let Some(k) = load_aead() {
        c.set_aead_key(k);
    }
    c
}

#[tokio::main(flavor = "multi_thread")]
async fn main() {
    let a: Vec<String> = std::env::args().collect();
    let server = a.get(1).cloned().unwrap_or_else(|| "http://127.0.0.1:50071".into());
    let model = a.get(2).cloned().unwrap_or_else(|| "full_alignment_ont_hac_v500".into());
    let batch: usize = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(64);
    let conc: usize = a.get(4).and_then(|s| s.parse().ok()).unwrap_or(64);
    let secs: u64 = a.get(5).and_then(|s| s.parse().ok()).unwrap_or(30);

    let mode = if tls_opts().is_some() {
        if std::env::var("FA_TLS_CERT").is_ok() { "mtls" } else { "tls" }
    } else {
        "plaintext"
    };
    let aead = if load_aead().is_some() { "+aead" } else { "" };

    let shape = vec![batch as i64, 89, 33, 8];
    let data = vec![1i8; batch * 89 * 33 * 8];
    let samples = Arc::new(AtomicU64::new(0));
    let reqs = Arc::new(AtomicU64::new(0));

    let t_connect = Instant::now();
    let mut probe = make_client(&server).await;
    let first_rt = {
        let t = Instant::now();
        let req = make_request_i8("probe", &model, shape.clone(), &data);
        probe.predict(req).await.expect("probe predict");
        t.elapsed().as_secs_f64() * 1000.0
    };
    let connect_plus_first = t_connect.elapsed().as_secs_f64() * 1000.0;

    let start = Instant::now();
    let deadline = start + Duration::from_secs(secs);
    let mut tasks = Vec::new();
    for _ in 0..conc {
        let (server, model, shape, data) = (server.clone(), model.clone(), shape.clone(), data.clone());
        let (samples, reqs) = (samples.clone(), reqs.clone());
        tasks.push(tokio::spawn(async move {
            let mut client = make_client(&server).await;
            while Instant::now() < deadline {
                let req = make_request_i8("load", &model, shape.clone(), &data);
                if client.predict(req).await.is_ok() {
                    samples.fetch_add(batch as u64, Ordering::Relaxed);
                    reqs.fetch_add(1, Ordering::Relaxed);
                }
            }
        }));
    }
    for t in tasks {
        let _ = t.await;
    }

    let el = start.elapsed().as_secs_f64();
    let s = samples.load(Ordering::Relaxed);
    let r = reqs.load(Ordering::Relaxed);
    let rps = r as f64 / el;
    let mean_lat_ms = if rps > 0.0 { conc as f64 / rps * 1000.0 } else { 0.0 };
    println!(
        "FA_LOAD mode={}{} model={} batch={} conc={} secs={:.1} | requests={} ({:.0}/s) | samples={} ({:.0} samples/s) | mean_lat={:.3}ms | first_rt={:.3}ms connect+first={:.3}ms",
        mode, aead, model, batch, conc, el, r, rps, s, s as f64 / el, mean_lat_ms, first_rt, connect_plus_first
    );
}
