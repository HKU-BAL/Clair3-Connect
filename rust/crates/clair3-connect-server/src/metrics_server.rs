use std::fmt::Write as _;
use std::net::SocketAddr;

use tokio::io::{AsyncReadExt, AsyncWriteExt};
use tokio::net::{TcpListener, TcpStream};

use crate::service::PredictService;

#[derive(Debug, Clone)]
pub struct ModelMetrics {

    pub model_id: String,

    pub model_sha256: String,

    pub device: String,

    pub batches: u64,

    pub items: u64,

    pub compute_ns: u64,

    pub gate_max_inflight: usize,

    pub gate_available_permits: usize,
}

#[derive(Debug, Clone)]
pub struct MetricsSnapshot {

    pub gate_max_inflight: usize,

    pub gate_available_permits: usize,

    pub per_model: Vec<ModelMetrics>,
}

impl MetricsSnapshot {

    pub fn gate_inflight(&self) -> usize {
        self.gate_max_inflight
            .saturating_sub(self.gate_available_permits)
    }
}

pub fn render_prometheus(s: &MetricsSnapshot) -> String {
    let mut out = String::with_capacity(1024);

    let _ = writeln!(
        out,
        "# HELP clair3_connect_gate_max_inflight Configured maximum concurrent in-flight Predict RPCs."
    );
    let _ = writeln!(out, "# TYPE clair3_connect_gate_max_inflight gauge");
    let _ = writeln!(out, "clair3_connect_gate_max_inflight {}", s.gate_max_inflight);

    let _ = writeln!(
        out,
        "# HELP clair3_connect_gate_inflight Currently in-flight Predict RPCs (max_inflight - available_permits)."
    );
    let _ = writeln!(out, "# TYPE clair3_connect_gate_inflight gauge");
    let _ = writeln!(out, "clair3_connect_gate_inflight {}", s.gate_inflight());

    let _ = writeln!(
        out,
        "# HELP clair3_connect_gate_available_permits Free permits on the backpressure semaphore."
    );
    let _ = writeln!(out, "# TYPE clair3_connect_gate_available_permits gauge");
    let _ = writeln!(
        out,
        "clair3_connect_gate_available_permits {}",
        s.gate_available_permits
    );

    let _ = writeln!(
        out,
        "# HELP clair3_connect_batches_total Cumulative count of batches processed by the dynamic batcher, per model."
    );
    let _ = writeln!(out, "# TYPE clair3_connect_batches_total counter");
    for m in &s.per_model {
        let _ = writeln!(
            out,
            "clair3_connect_batches_total{{model_id=\"{}\",model_sha256=\"{}\",device=\"{}\"}} {}",
            escape_label(&m.model_id),
            escape_label(&m.model_sha256),
            escape_label(&m.device),
            m.batches
        );
    }

    let _ = writeln!(
        out,
        "# HELP clair3_connect_items_total Cumulative count of items (rows) processed by the dynamic batcher, per model and device."
    );
    let _ = writeln!(out, "# TYPE clair3_connect_items_total counter");
    for m in &s.per_model {
        let _ = writeln!(
            out,
            "clair3_connect_items_total{{model_id=\"{}\",model_sha256=\"{}\",device=\"{}\"}} {}",
            escape_label(&m.model_id),
            escape_label(&m.model_sha256),
            escape_label(&m.device),
            m.items
        );
    }

    let _ = writeln!(
        out,
        "# HELP clair3_connect_compute_ns_total Cumulative nanoseconds of server-side compute (cat+forward+split) per model and device."
    );
    let _ = writeln!(out, "# TYPE clair3_connect_compute_ns_total counter");
    for m in &s.per_model {
        let _ = writeln!(
            out,
            "clair3_connect_compute_ns_total{{model_id=\"{}\",model_sha256=\"{}\",device=\"{}\"}} {}",
            escape_label(&m.model_id),
            escape_label(&m.model_sha256),
            escape_label(&m.device),
            m.compute_ns
        );
    }

    let _ = writeln!(
        out,
        "# HELP clair3_connect_per_device_gate_max_inflight Per-device backpressure cap (max_inflight)."
    );
    let _ = writeln!(out, "# TYPE clair3_connect_per_device_gate_max_inflight gauge");
    for m in &s.per_model {
        let _ = writeln!(
            out,
            "clair3_connect_per_device_gate_max_inflight{{model_id=\"{}\",device=\"{}\"}} {}",
            escape_label(&m.model_id),
            escape_label(&m.device),
            m.gate_max_inflight
        );
    }
    let _ = writeln!(
        out,
        "# HELP clair3_connect_per_device_gate_available_permits Free permits on each device's backpressure semaphore."
    );
    let _ = writeln!(out, "# TYPE clair3_connect_per_device_gate_available_permits gauge");
    for m in &s.per_model {
        let _ = writeln!(
            out,
            "clair3_connect_per_device_gate_available_permits{{model_id=\"{}\",device=\"{}\"}} {}",
            escape_label(&m.model_id),
            escape_label(&m.device),
            m.gate_available_permits
        );
    }

    out
}

fn escape_label(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            other => out.push(other),
        }
    }
    out
}

pub async fn spawn(
    addr: SocketAddr,
    svc: PredictService,
) -> std::io::Result<tokio::task::JoinHandle<()>> {
    let listener = TcpListener::bind(addr).await?;
    tracing::info!(
        bind = %addr,
        "metrics server listening on /metrics and /healthz"
    );
    let handle = tokio::spawn(async move {
        loop {
            let (stream, _peer) = match listener.accept().await {
                Ok(p) => p,
                Err(e) => {
                    tracing::warn!("metrics accept: {}", e);
                    continue;
                }
            };
            let svc = svc.clone();
            tokio::spawn(async move {
                if let Err(e) = handle_connection(stream, svc).await {
                    tracing::debug!("metrics conn: {}", e);
                }
            });
        }
    });
    Ok(handle)
}

const READ_TIMEOUT: std::time::Duration = std::time::Duration::from_secs(2);

async fn handle_connection(mut stream: TcpStream, svc: PredictService) -> std::io::Result<()> {
    let mut buf = [0u8; 4096];
    let mut total = 0usize;

    let read_fut = async {
        loop {
            let n = stream.read(&mut buf[total..]).await?;
            if n == 0 {
                break Ok::<_, std::io::Error>(());
            }
            total += n;
            if buf[..total].windows(4).any(|w| w == b"\r\n\r\n") {
                break Ok(());
            }
            if total == buf.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "headers exceed 4 KiB",
                ));
            }
        }
    };
    match tokio::time::timeout(READ_TIMEOUT, read_fut).await {
        Ok(Ok(())) => {}
        Ok(Err(e)) if e.kind() == std::io::ErrorKind::InvalidData => {
            return write_response(&mut stream, 413, "Payload Too Large", "").await;
        }
        Ok(Err(e)) => return Err(e),
        Err(_) => {
            return write_response(&mut stream, 408, "Request Timeout", "").await;
        }
    }
    let head = std::str::from_utf8(&buf[..total]).unwrap_or("");
    let first = head.lines().next().unwrap_or("");
    let mut parts = first.split_whitespace();
    let method = parts.next().unwrap_or("");
    let path = parts.next().unwrap_or("");
    if method != "GET" {
        return write_response(&mut stream, 405, "Method Not Allowed", "").await;
    }
    match path {
        "/metrics" => {
            let body = render_prometheus(&svc.metrics_snapshot());
            write_full_response(
                &mut stream,
                200,
                "OK",
                "text/plain; version=0.0.4; charset=utf-8",
                &body,
            )
            .await
        }
        "/healthz" => {
            write_full_response(&mut stream, 200, "OK", "text/plain; charset=utf-8", "ok\n")
                .await
        }
        _ => write_response(&mut stream, 404, "Not Found", "").await,
    }
}

async fn write_response(
    stream: &mut TcpStream,
    code: u16,
    reason: &str,
    body: &str,
) -> std::io::Result<()> {
    write_full_response(stream, code, reason, "text/plain; charset=utf-8", body).await
}

async fn write_full_response(
    stream: &mut TcpStream,
    code: u16,
    reason: &str,
    content_type: &str,
    body: &str,
) -> std::io::Result<()> {
    let resp = format!(
        "HTTP/1.1 {} {}\r\nContent-Type: {}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{}",
        code,
        reason,
        content_type,
        body.len(),
        body
    );
    stream.write_all(resp.as_bytes()).await?;
    stream.shutdown().await
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn sample() -> MetricsSnapshot {
        MetricsSnapshot {
            gate_max_inflight: 8,
            gate_available_permits: 5,
            per_model: vec![
                ModelMetrics {
                    model_id: "pileup_ont_hac_v500".to_string(),
                    model_sha256: "deadbeef".to_string(),
                    device: "cuda:0".to_string(),
                    batches: 42,
                    items: 1337,
                    compute_ns: 123_456,
                    gate_max_inflight: 8,
                    gate_available_permits: 5,
                },
                ModelMetrics {
                    model_id: "fa_ont_hac".to_string(),
                    model_sha256: "cafef00d".to_string(),
                    device: "cuda:0".to_string(),
                    batches: 7,
                    items: 88,
                    compute_ns: 7_000,
                    gate_max_inflight: 8,
                    gate_available_permits: 8,
                },
            ],
        }
    }

    #[test]
    fn snapshot_inflight_calculation() {
        let s = sample();
        assert_eq!(s.gate_inflight(), 3, "8 max - 5 free = 3 inflight");
    }

    #[test]
    fn snapshot_inflight_saturates() {

        let s = MetricsSnapshot {
            gate_max_inflight: 4,
            gate_available_permits: 10,
            per_model: vec![],
        };
        assert_eq!(s.gate_inflight(), 0);
    }

    #[test]
    fn render_contains_expected_series() {
        let body = render_prometheus(&sample());

        assert!(body.contains("# TYPE clair3_connect_gate_max_inflight gauge"));
        assert!(body.contains("# TYPE clair3_connect_gate_inflight gauge"));
        assert!(body.contains("# TYPE clair3_connect_gate_available_permits gauge"));
        assert!(body.contains("# TYPE clair3_connect_batches_total counter"));
        assert!(body.contains("# TYPE clair3_connect_items_total counter"));
        assert!(body.contains("# TYPE clair3_connect_compute_ns_total counter"));

        assert!(body.contains("clair3_connect_gate_max_inflight 8"));
        assert!(body.contains("clair3_connect_gate_available_permits 5"));
        assert!(body.contains("clair3_connect_gate_inflight 3"));

        assert!(body.contains(
            "clair3_connect_batches_total{model_id=\"pileup_ont_hac_v500\",model_sha256=\"deadbeef\",device=\"cuda:0\"} 42"
        ));
        assert!(body.contains(
            "clair3_connect_items_total{model_id=\"pileup_ont_hac_v500\",model_sha256=\"deadbeef\",device=\"cuda:0\"} 1337"
        ));
        assert!(body.contains(
            "clair3_connect_compute_ns_total{model_id=\"pileup_ont_hac_v500\",model_sha256=\"deadbeef\",device=\"cuda:0\"} 123456"
        ));

        assert!(body.contains(
            "clair3_connect_per_device_gate_max_inflight{model_id=\"pileup_ont_hac_v500\",device=\"cuda:0\"} 8"
        ));
        assert!(body.contains(
            "clair3_connect_per_device_gate_available_permits{model_id=\"pileup_ont_hac_v500\",device=\"cuda:0\"} 5"
        ));
    }

    #[test]
    fn escape_label_escapes_problematic_chars() {
        assert_eq!(escape_label("a\\b"), "a\\\\b");
        assert_eq!(escape_label("a\"b"), "a\\\"b");
        assert_eq!(escape_label("a\nb"), "a\\nb");
        assert_eq!(escape_label("ok"), "ok");
    }

    #[test]
    fn no_duplicate_series_per_metric() {
        let body = render_prometheus(&sample());
        let mut seen: HashMap<String, usize> = HashMap::new();
        for line in body.lines() {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }

            let (key, _value) = line.rsplit_once(' ').unwrap_or((line, ""));
            *seen.entry(key.to_string()).or_insert(0) += 1;
        }
        for (k, v) in &seen {
            assert_eq!(*v, 1, "duplicate series for {k}");
        }
    }

    #[test]
    #[should_panic(expected = "duplicate series")]
    fn dedup_test_catches_actual_duplicates() {
        let body = "# TYPE x counter\nx{a=\"b\"} 1\nx{a=\"b\"} 2\n";
        let mut seen: HashMap<String, usize> = HashMap::new();
        for line in body.lines() {
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let (key, _value) = line.rsplit_once(' ').unwrap_or((line, ""));
            *seen.entry(key.to_string()).or_insert(0) += 1;
        }
        for (k, v) in &seen {
            assert_eq!(*v, 1, "duplicate series for {k}");
        }
    }
}
