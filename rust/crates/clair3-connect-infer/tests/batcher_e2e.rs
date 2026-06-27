use clair3_connect_infer::batcher::{BatchConfig, BatchedWorker};
use clair3_connect_infer::model_registry::Registry;
use clair3_connect_infer::torchscript::LoadedModel;
use clair3_connect_infer::worker::Worker;
use std::path::PathBuf;
use std::sync::atomic::Ordering;
use std::time::Duration;
use tch::Tensor;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

#[tokio::test(flavor = "multi_thread", worker_threads = 4)]
async fn concurrent_requests_get_batched() {
    let yaml = workspace_root().join("models/model_registry.example.yaml");
    if !yaml.exists() {
        eprintln!("Skipping — registry yaml missing");
        return;
    }
    let reg = Registry::load(&yaml).expect("registry");
    let yaml_dir = yaml.parent().unwrap().to_path_buf();
    let entry = reg
        .models
        .iter()
        .find(|e| e.model_id == "pileup_ont_hac_v500")
        .expect("pileup_ont_hac_v500 in registry");
    let ts_path = match reg.verify_torchscript(entry, &yaml_dir) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("Skipping — TorchScript verify failed: {e}");
            return;
        }
    };
    let model = match LoadedModel::load(&ts_path, entry.clone(), tch::Device::Cpu) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("Skipping — model load failed: {e}");
            return;
        }
    };
    let worker = Worker::new(model);

    let cfg = BatchConfig {
        max_batch_size: 64,
        flush_timeout: Duration::from_millis(32),
        queue_depth: 256,
    };
    let bw = BatchedWorker::spawn(worker, cfg);

    let n: usize = 16;
    let mut handles = Vec::with_capacity(n);
    for _ in 0..n {
        let bw = bw.clone();
        handles.push(tokio::spawn(async move {
            let x = Tensor::zeros([1, 33, 18], (tch::Kind::Int, tch::Device::Cpu));
            bw.predict(x).await.expect("predict")
        }));
    }
    let mut outputs = Vec::with_capacity(n);
    for h in handles {
        outputs.push(h.await.unwrap());
    }

    for y in &outputs {
        let size = y.size();
        assert_eq!(size.len(), 2, "y shape rank");
        assert_eq!(size[0], 1, "y batch dim");
        assert_eq!(size[1], 24, "y output width");
    }

    let m = bw.metrics();
    let items = m.items.load(Ordering::Relaxed);
    let batches = m.batches.load(Ordering::Relaxed);
    eprintln!("batcher: items={items}, batches={batches}");
    assert_eq!(items, n as u64, "all items processed");
    assert!(
        batches < n as u64,
        "batcher should consolidate concurrent requests: batches={batches} >= n={n}"
    );

    assert!(
        batches <= (n as u64 / 2 + 1),
        "expected at most ~n/2 batches for 16 concurrent requests with a 32 ms window; got {batches}"
    );
}
