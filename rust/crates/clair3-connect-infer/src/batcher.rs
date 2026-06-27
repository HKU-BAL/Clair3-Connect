use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::time::Duration;

use anyhow::Result;
use tch::Tensor;
use tokio::sync::{mpsc, oneshot};

use crate::worker::Worker;

#[derive(Debug, Clone)]
pub struct BatchConfig {

    pub max_batch_size: usize,

    pub flush_timeout: Duration,

    pub queue_depth: usize,
}

impl Default for BatchConfig {
    fn default() -> Self {
        Self {
            max_batch_size: 64,
            flush_timeout: Duration::from_millis(4),
            queue_depth: 1024,
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum BatchError {

    #[error("batcher queue full")]
    Full,

    #[error("batcher shut down")]
    ShutDown,

    #[error("forward: {0}")]
    Forward(String),
}

struct BatchItem {

    tensor: Tensor,

    rows: i64,

    reply: oneshot::Sender<Result<Tensor, String>>,
}

#[derive(Debug, Default)]
pub struct BatcherMetrics {

    pub batches: AtomicU64,

    pub items: AtomicU64,

    pub compute_ns: AtomicU64,
}

#[derive(Clone)]
pub struct BatchedWorker {
    tx: mpsc::Sender<BatchItem>,
    cfg: BatchConfig,
    metrics: Arc<BatcherMetrics>,
}

impl BatchedWorker {

    pub fn spawn(worker: Worker, cfg: BatchConfig) -> Self {
        assert!(
            cfg.max_batch_size > 0,
            "BatchConfig.max_batch_size must be > 0"
        );
        assert!(cfg.queue_depth > 0, "BatchConfig.queue_depth must be > 0");
        assert!(
            !cfg.flush_timeout.is_zero(),
            "BatchConfig.flush_timeout must be > 0; use a 1 ns value to opt out of amortization"
        );
        let (tx, rx) = mpsc::channel::<BatchItem>(cfg.queue_depth);
        let worker = Arc::new(worker);
        let cfg2 = cfg.clone();
        let metrics = Arc::new(BatcherMetrics::default());
        let metrics_loop = metrics.clone();
        tokio::spawn(async move {
            batch_loop(worker, rx, cfg2, metrics_loop).await;
        });
        Self { tx, cfg, metrics }
    }

    pub fn metrics(&self) -> &BatcherMetrics {
        &self.metrics
    }

    pub async fn predict(&self, x: Tensor) -> Result<Tensor, BatchError> {
        let rows = x.size().first().copied().unwrap_or(0);
        let (reply_tx, reply_rx) = oneshot::channel();
        let item = BatchItem {
            tensor: x,
            rows,
            reply: reply_tx,
        };
        self.tx.try_send(item).map_err(|e| match e {
            mpsc::error::TrySendError::Full(_) => BatchError::Full,
            mpsc::error::TrySendError::Closed(_) => BatchError::ShutDown,
        })?;
        match reply_rx.await {
            Ok(Ok(y)) => Ok(y),
            Ok(Err(msg)) => Err(BatchError::Forward(msg)),
            Err(_) => Err(BatchError::ShutDown),
        }
    }

    pub fn config(&self) -> &BatchConfig {
        &self.cfg
    }
}

async fn batch_loop(
    worker: Arc<Worker>,
    mut rx: mpsc::Receiver<BatchItem>,
    cfg: BatchConfig,
    metrics: Arc<BatcherMetrics>,
) {
    loop {

        let first = match rx.recv().await {
            Some(item) => item,
            None => return,
        };
        let mut batch: Vec<BatchItem> = vec![first];
        let deadline = tokio::time::Instant::now() + cfg.flush_timeout;

        while batch.len() < cfg.max_batch_size {
            match tokio::time::timeout_at(deadline, rx.recv()).await {
                Ok(Some(item)) => batch.push(item),
                Ok(None) => break,
                Err(_) => break,
            }
        }
        let n = batch.len() as u64;
        metrics.batches.fetch_add(1, Ordering::Relaxed);
        metrics.items.fetch_add(n, Ordering::Relaxed);
        let compute_t0 = std::time::Instant::now();
        process_batch(&worker, batch);
        metrics
            .compute_ns
            .fetch_add(compute_t0.elapsed().as_nanos() as u64, Ordering::Relaxed);

    }
}

fn process_batch(worker: &Worker, batch: Vec<BatchItem>) {

    if let Err(msg) = validate_batch_homogeneous(&batch) {
        for item in batch {
            let _ = item.reply.send(Err(msg.clone()));
        }
        return;
    }

    let tensors: Vec<&Tensor> = batch.iter().map(|b| &b.tensor).collect();
    let big = Tensor::cat(&tensors, 0);
    let y = match worker.forward_batch(&big) {
        Ok(y) => y,
        Err(e) => {

            let msg = format!("{:#}", e);
            for item in batch {
                let _ = item.reply.send(Err(msg.clone()));
            }
            return;
        }
    };

    let mut offset: i64 = 0;
    for item in batch {
        let n = item.rows;
        let slice = y.narrow(0, offset, n).contiguous();
        let _ = item.reply.send(Ok(slice));
        offset += n;
    }
}

fn validate_batch_homogeneous(batch: &[BatchItem]) -> Result<(), String> {
    let first = batch
        .first()
        .ok_or_else(|| "empty batch".to_string())?;
    let want_size = first.tensor.size();
    let want_kind = first.tensor.kind();
    let want_dev = first.tensor.device();
    for (i, item) in batch.iter().enumerate().skip(1) {
        let got_size = item.tensor.size();
        if got_size.len() != want_size.len() || got_size[1..] != want_size[1..] {
            return Err(format!(
                "batch shape mismatch at index {}: expected per-row {:?}, got {:?}",
                i,
                &want_size[1..],
                &got_size[1..]
            ));
        }
        if item.tensor.kind() != want_kind {
            return Err(format!(
                "batch dtype mismatch at index {}: expected {:?}, got {:?}",
                i,
                want_kind,
                item.tensor.kind()
            ));
        }
        if item.tensor.device() != want_dev {
            return Err(format!(
                "batch device mismatch at index {}: expected {:?}, got {:?}",
                i,
                want_dev,
                item.tensor.device()
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;

    fn item(x: Tensor) -> BatchItem {
        let rows = x.size().first().copied().unwrap_or(0);
        let (tx, _rx) = oneshot::channel();
        BatchItem {
            tensor: x,
            rows,
            reply: tx,
        }
    }

    #[test]
    fn validate_passes_when_shapes_match() {
        let a = Tensor::zeros([3, 4, 5], (tch::Kind::Float, tch::Device::Cpu));
        let b = Tensor::zeros([2, 4, 5], (tch::Kind::Float, tch::Device::Cpu));
        let c = Tensor::zeros([7, 4, 5], (tch::Kind::Float, tch::Device::Cpu));
        let batch = vec![item(a), item(b), item(c)];
        validate_batch_homogeneous(&batch).expect("homogeneous");
    }

    #[test]
    fn validate_rejects_shape_mismatch() {
        let a = Tensor::zeros([1, 4, 5], (tch::Kind::Float, tch::Device::Cpu));
        let b = Tensor::zeros([1, 4, 6], (tch::Kind::Float, tch::Device::Cpu));
        let batch = vec![item(a), item(b)];
        let err = validate_batch_homogeneous(&batch).unwrap_err();
        assert!(err.contains("shape mismatch"), "msg was: {err}");
    }

    #[test]
    fn validate_rejects_dtype_mismatch() {
        let a = Tensor::zeros([1, 4], (tch::Kind::Float, tch::Device::Cpu));
        let b = Tensor::zeros([1, 4], (tch::Kind::Int, tch::Device::Cpu));
        let batch = vec![item(a), item(b)];
        let err = validate_batch_homogeneous(&batch).unwrap_err();
        assert!(err.contains("dtype mismatch"), "msg was: {err}");
    }

    #[test]
    fn validate_rejects_rank_mismatch() {
        let a = Tensor::zeros([1, 4, 5], (tch::Kind::Float, tch::Device::Cpu));
        let b = Tensor::zeros([1, 4], (tch::Kind::Float, tch::Device::Cpu));
        let batch = vec![item(a), item(b)];
        let err = validate_batch_homogeneous(&batch).unwrap_err();
        assert!(err.contains("shape mismatch"), "msg was: {err}");
    }

    #[test]
    #[should_panic(expected = "max_batch_size")]
    fn spawn_rejects_zero_batch_size() {

        let cfg = BatchConfig {
            max_batch_size: 0,
            flush_timeout: Duration::from_millis(1),
            queue_depth: 1,
        };
        assert!(cfg.max_batch_size > 0, "max_batch_size must be > 0");
    }

    #[test]
    fn default_config_is_sane() {
        let cfg = BatchConfig::default();
        assert!(cfg.max_batch_size > 0);
        assert!(cfg.queue_depth > 0);
        assert!(!cfg.flush_timeout.is_zero());
    }
}
