use std::path::PathBuf;

use tch::{Device, Kind, Tensor};

use clair3_connect_infer::model_registry::Registry;
use clair3_connect_infer::torchscript::LoadedModel;
use clair3_connect_infer::worker::Worker;

const MODEL_ID: &str = "full_alignment_ont_hac_v500";

fn registry_yaml() -> PathBuf {
    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    crate_dir
        .join("..")
        .join("..")
        .join("models")
        .join("model_registry.example.yaml")
}

fn load_fa(device: Device) -> Option<Worker> {
    let yaml = registry_yaml();
    if !yaml.exists() {
        eprintln!("[skip] registry yaml not at {}", yaml.display());
        return None;
    }
    let reg = Registry::load(&yaml).expect("parse registry");
    let entry = reg.get(MODEL_ID).expect("FA model registered");
    let ts = match reg.verify_torchscript(entry, yaml.parent().unwrap()) {
        Ok(p) => p,
        Err(e) => {
            if e.to_string().contains("No such file") {
                eprintln!("[skip] {}", e);
                return None;
            }
            panic!("registry verify failed: {e}");
        }
    };
    let model = LoadedModel::load(&ts, entry.clone(), device).expect("FA CModule load");
    Some(Worker::new(model))
}

fn softmax_heads_sum_to_one(y: &[f32], n: usize, c: usize) -> Option<(usize, [f32; 4])> {

    assert_eq!(c, 90, "FA output should be 90 wide");
    let heads = [(0_usize, 21_usize), (21, 24), (24, 57), (57, 90)];
    for row in 0..n {
        let r = row * c;
        let sums: [f32; 4] = heads.map(|(a, b)| y[r + a..r + b].iter().sum::<f32>());
        for &s in &sums {
            if (s - 1.0).abs() > 1e-3 {
                return Some((row, sums));
            }
        }
    }
    None
}

#[test]
fn fa_smoke_load_and_forward() {
    let device = Device::Cpu;
    let Some(worker) = load_fa(device) else {
        return;
    };

    let entry = &worker.model.entry;
    assert_eq!(entry.input.shape, vec![-1, 89, 33, 8]);
    assert_eq!(entry.input.dtype, "int8");
    assert_eq!(entry.output.shape, vec![-1, 90]);

    for &batch in &[1, 2, 4, 8] {
        let x = Tensor::zeros([batch as i64, 89, 33, 8], (Kind::Int8, device));
        let y = worker.forward_batch(&x).expect("FA forward");
        let yshape = y.size();
        assert_eq!(
            yshape,
            vec![batch as i64, 90],
            "FA batch {} output shape",
            batch
        );
        assert_eq!(y.kind(), Kind::Float, "FA output should be float32");

        let y_flat = y
            .to_kind(Kind::Float)
            .contiguous()
            .reshape([batch as i64 * 90]);
        let v: Vec<f32> = Vec::<f32>::try_from(&y_flat).expect("Y -> Vec");
        if let Some((row, sums)) = softmax_heads_sum_to_one(&v, batch, 90) {
            panic!(
                "FA batch {} row {} head sums {:?} not all ≈ 1.0",
                batch, row, sums
            );
        }
        eprintln!(
            "[fa-smoke] batch={} shape=[{},90] head sums OK",
            batch, batch
        );
    }
}
