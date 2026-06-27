use std::fs::File;
use std::path::{Path, PathBuf};

use ndarray::{Array2, Array3};
use ndarray_npy::ReadNpyExt;
use tch::{Device, Kind, Tensor};

use clair3_connect_infer::device::Device as Clair3Device;
use clair3_connect_infer::model_registry::Registry;
use clair3_connect_infer::torchscript::LoadedModel;
use clair3_connect_infer::worker::Worker;

const MODEL_ID: &str = "pileup_ont_hac_v500";

fn golden_dir() -> Option<PathBuf> {
    std::env::var_os("CLAIR3_CONNECT_GOLDEN_DIR").map(PathBuf::from)
}

fn registry_yaml() -> PathBuf {
    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    crate_dir
        .join("..")
        .join("..")
        .join("models")
        .join("model_registry.example.yaml")
}

fn load_pileup(device: tch::Device) -> Option<Worker> {
    let yaml = registry_yaml();
    if !yaml.exists() {
        eprintln!("[skip] registry yaml not at {}", yaml.display());
        return None;
    }
    let reg = Registry::load(&yaml).expect("parse registry");
    let entry = reg.get(MODEL_ID).expect("pileup model registered");
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
    let model = LoadedModel::load(&ts, entry.clone(), device).expect("CModule load");
    Some(Worker::new(model))
}

fn load_phase0_pileup(dir: &Path, y_name: &str) -> Option<(Array3<i32>, Array2<f32>)> {
    let x_path = dir.join("pileup_post_scale_X_int32.npy");
    let y_path = dir.join(y_name);
    if !x_path.exists() || !y_path.exists() {
        eprintln!(
            "[skip] missing fixture: {} or {}",
            x_path.display(),
            y_path.display()
        );
        return None;
    }
    let x = Array3::<i32>::read_npy(File::open(&x_path).unwrap()).expect("read X");
    let y = Array2::<f32>::read_npy(File::open(&y_path).unwrap()).expect("read Y");
    Some((x, y))
}

fn tensor_from_array3_i32(arr: &Array3<i32>, device: Device) -> Tensor {
    let s = arr.shape();
    let shape: Vec<i64> = s.iter().map(|&d| d as i64).collect();
    let slice = arr.as_slice().expect("contiguous");
    Tensor::from_slice(slice)
        .view(shape.as_slice())
        .to_device(device)
}

fn diff_stats(y_rust: &Tensor, y_py: &Array2<f32>) -> (f32, f32) {
    let y_rust_cpu = y_rust
        .to_device(Device::Cpu)
        .to_kind(Kind::Float)
        .contiguous();
    let n = y_py.shape()[0];
    let c = y_py.shape()[1];
    let total = (n * c) as i64;
    let flat = y_rust_cpu.reshape([total]);
    let v: Vec<f32> = Vec::<f32>::try_from(&flat).expect("Y rust -> Vec<f32>");
    assert_eq!(
        v.len(),
        n * c,
        "Y rust size {} != Y python size {}",
        v.len(),
        n * c
    );
    let py = y_py.as_slice().expect("contig py");
    let mut max_abs = 0.0_f32;
    let mut argmax_eq_rows = 0usize;
    for row in 0..n {
        let r = row * c;
        for k in 0..c {
            let d = (v[r + k] - py[r + k]).abs();
            if d > max_abs {
                max_abs = d;
            }
        }

        let am_rust_gt21 = v[r..r + 21]
            .iter()
            .enumerate()
            .fold((0usize, f32::NEG_INFINITY), |(am, mx), (i, &x)| {
                if x > mx {
                    (i, x)
                } else {
                    (am, mx)
                }
            })
            .0;
        let am_py_gt21 = py[r..r + 21]
            .iter()
            .enumerate()
            .fold((0usize, f32::NEG_INFINITY), |(am, mx), (i, &x)| {
                if x > mx {
                    (i, x)
                } else {
                    (am, mx)
                }
            })
            .0;
        let am_rust_gt = v[r + 21..r + 24]
            .iter()
            .enumerate()
            .fold((0usize, f32::NEG_INFINITY), |(am, mx), (i, &x)| {
                if x > mx {
                    (i, x)
                } else {
                    (am, mx)
                }
            })
            .0;
        let am_py_gt = py[r + 21..r + 24]
            .iter()
            .enumerate()
            .fold((0usize, f32::NEG_INFINITY), |(am, mx), (i, &x)| {
                if x > mx {
                    (i, x)
                } else {
                    (am, mx)
                }
            })
            .0;
        if am_rust_gt21 == am_py_gt21 && am_rust_gt == am_py_gt {
            argmax_eq_rows += 1;
        }
    }
    let frac = argmax_eq_rows as f32 / n as f32;
    (max_abs, frac)
}

#[test]
fn pileup_y_parity_cpu() {
    let Some(dir) = golden_dir() else {
        eprintln!("[skip] CLAIR3_CONNECT_GOLDEN_DIR not set");
        return;
    };
    let Some((x_np, y_py)) = load_phase0_pileup(&dir, "Y_pileup_cpu.npy") else {
        return;
    };

    let device = Clair3Device::Cpu.as_tch();
    let Some(worker) = load_pileup(device) else {
        return;
    };

    let x = tensor_from_array3_i32(&x_np, device);
    let y_rust = worker.forward_batch(&x).expect("forward");

    let (max_abs, argmax_match) = diff_stats(&y_rust, &y_py);
    eprintln!(
        "[parity-cpu] N={} max_abs_diff={:.3e}  argmax_match={:.4}",
        x_np.shape()[0],
        max_abs,
        argmax_match
    );

    let tol = 1e-5_f32;
    assert!(
        max_abs < tol,
        "CPU max_abs_diff {} exceeds tolerance {}",
        max_abs,
        tol
    );
    assert!(
        argmax_match >= 1.0 - 1e-9,
        "CPU argmax match {} < 1.0 (expected exact)",
        argmax_match
    );
}

#[test]
fn pileup_y_parity_cuda() {
    let Some(dir) = golden_dir() else {
        eprintln!("[skip] CLAIR3_CONNECT_GOLDEN_DIR not set");
        return;
    };
    if !tch::Cuda::is_available() {
        eprintln!("[skip] CUDA not available");
        return;
    }
    let Some((x_np, y_py)) = load_phase0_pileup(&dir, "Y_pileup_cuda.npy") else {
        return;
    };

    let device = Clair3Device::Cuda(0).as_tch();
    let Some(worker) = load_pileup(device) else {
        return;
    };

    let x = tensor_from_array3_i32(&x_np, device);
    let y_rust = worker.forward_batch(&x).expect("forward");

    let (max_abs, argmax_match) = diff_stats(&y_rust, &y_py);
    eprintln!(
        "[parity-cuda] N={} max_abs_diff={:.3e}  argmax_match={:.4}",
        x_np.shape()[0],
        max_abs,
        argmax_match
    );

    let tol = 1e-4_f32;
    assert!(
        max_abs < tol,
        "CUDA max_abs_diff {} exceeds tolerance {}",
        max_abs,
        tol
    );
    assert!(
        argmax_match >= 1.0 - 1e-9,
        "CUDA argmax match {} < 1.0 (expected exact)",
        argmax_match
    );
}
