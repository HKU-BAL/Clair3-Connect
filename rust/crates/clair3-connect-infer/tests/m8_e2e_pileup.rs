use std::fs::File;
use std::path::PathBuf;

use ndarray::Array2;
use ndarray_npy::ReadNpyExt;
use tch::{Device, Kind, Tensor};

use clair3_connect_core::hts::region::Region;
use clair3_connect_core::pileup::candidates::window_candidates;
use clair3_connect_core::pileup::kernel::{calculate_pileup, PileupParams};

use clair3_connect_infer::model_registry::Registry;
use clair3_connect_infer::torchscript::LoadedModel;
use clair3_connect_infer::worker::Worker;

const FIXTURE_REGION: &str = "chr20:99966-300033";
const PILEUP_MODEL_ID: &str = "pileup_ont_hac_v500";

fn golden_dir() -> Option<PathBuf> {
    std::env::var_os("CLAIR3_CONNECT_GOLDEN_DIR").map(PathBuf::from)
}

fn fixture_bam() -> PathBuf {
    workspace_root()
        .join("test_data")
        .join("ont_chr20_demo")
        .join("HG003_chr20_demo.bam")
}

fn fixture_ref() -> PathBuf {
    workspace_root()
        .join("test_data")
        .join("ont_chr20_demo")
        .join("GRCh38_no_alt_chr20.fa")
}

fn workspace_root() -> PathBuf {
    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    crate_dir.join("..").join("..")
}

fn registry_yaml() -> PathBuf {
    workspace_root()
        .join("models")
        .join("model_registry.example.yaml")
}

#[test]
fn m8_rust_pipeline_y_parity() {
    let Some(dir) = golden_dir() else {
        eprintln!("[skip] CLAIR3_CONNECT_GOLDEN_DIR not set");
        return;
    };
    let y_path = dir.join("Y_pileup_cpu.npy");
    if !y_path.exists() {
        eprintln!("[skip] {} not present", y_path.display());
        return;
    }
    if !fixture_bam().exists() || !fixture_ref().exists() {
        eprintln!("[skip] missing fixture BAM/ref");
        return;
    }

    let region = Region::parse(FIXTURE_REGION).expect("parse region");
    let params = PileupParams::default();
    let res = calculate_pileup(&region, &fixture_bam(), &fixture_ref(), &params)
        .expect("calculate_pileup");
    let windowed = window_candidates(&res.chunks, &res.candidates);
    let n = windowed.len();
    eprintln!("[m8] Rust pipeline produced {n} windowed candidates");

    let max_depth: i64 = params.max_depth as i64;
    let scale_threshold = max_depth * 3 / 2;
    let mut n_scaled = 0;
    let mut flat: Vec<i32> = Vec::with_capacity(n * 33 * 18);
    for wc in &windowed {

        let depth: i64 = wc
            .alt_info
            .split('-')
            .next()
            .and_then(|s| s.parse().ok())
            .expect("alt_info has depth prefix");
        if depth > 0 && depth > scale_threshold {
            n_scaled += 1;
            let scale = depth as f64 / max_depth as f64;
            for &v in &wc.tensor {
                let scaled = (v as f64 / scale) as i32;
                flat.push(scaled);
            }
        } else {
            flat.extend_from_slice(&wc.tensor);
        }
    }
    eprintln!("[m8] {n_scaled} candidate(s) had depth scaling applied");
    assert_eq!(flat.len(), n * 33 * 18);

    let device = Device::Cpu;
    let x_tensor = Tensor::from_slice(&flat)
        .view([n as i64, 33, 18])
        .to_device(device);
    assert_eq!(x_tensor.kind(), Kind::Int, "X dtype should be int32");

    let yaml = registry_yaml();
    if !yaml.exists() {
        eprintln!("[skip] registry not at {}", yaml.display());
        return;
    }
    let reg = Registry::load(&yaml).expect("registry");
    let entry = reg.get(PILEUP_MODEL_ID).expect("pileup model registered");
    let ts_path = match reg.verify_torchscript(entry, yaml.parent().unwrap()) {
        Ok(p) => p,
        Err(e) => {
            if e.to_string().contains("No such file") {
                eprintln!("[skip] {}", e);
                return;
            }
            panic!("registry verify failed: {e}");
        }
    };
    let model = LoadedModel::load(&ts_path, entry.clone(), device).expect("load CModule");
    let worker = Worker::new(model);

    let y_rust = worker.forward_batch(&x_tensor).expect("forward");
    assert_eq!(y_rust.size(), vec![n as i64, 24], "Y rust shape mismatch");

    let y_py = Array2::<f32>::read_npy(File::open(&y_path).unwrap()).expect("read Y");
    assert_eq!(
        (y_py.shape()[0], y_py.shape()[1]),
        (n, 24),
        "candidate count mismatch: rust {} vs python {}",
        n,
        y_py.shape()[0]
    );

    let y_rust_cpu = y_rust.to_kind(Kind::Float).contiguous();
    let total = (n * 24) as i64;
    let y_rust_flat = y_rust_cpu.reshape([total]);
    let v: Vec<f32> = Vec::<f32>::try_from(&y_rust_flat).expect("Y rust → Vec");
    let py_slice = y_py.as_slice().expect("py contig");

    assert!(
        v.iter().all(|x| x.is_finite()),
        "Rust Y contains non-finite values (NaN/Inf)"
    );
    assert!(
        py_slice.iter().all(|x| x.is_finite()),
        "Python Y contains non-finite values (NaN/Inf)"
    );

    let mut max_abs = 0.0_f32;
    let mut argmax_eq = 0usize;
    for i in 0..n {
        let r = i * 24;
        for k in 0..24 {
            let d = (v[r + k] - py_slice[r + k]).abs();
            if d > max_abs {
                max_abs = d;
            }
        }

        let gt21_r = v[r..r + 21]
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
        let gt21_p = py_slice[r..r + 21]
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
        let gt_r = v[r + 21..r + 24]
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
        let gt_p = py_slice[r + 21..r + 24]
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
        if gt21_r == gt21_p && gt_r == gt_p {
            argmax_eq += 1;
        }
    }
    eprintln!(
        "[m8] N={n}, max_abs_diff={max_abs:.3e}, argmax_match={}/{n}",
        argmax_eq
    );

    let tol = 1e-5_f32;
    assert!(
        max_abs < tol,
        "M8 max_abs_diff {} exceeds tolerance {}",
        max_abs,
        tol
    );
    assert_eq!(
        argmax_eq,
        n,
        "M8 argmax mismatch on {} of {} rows",
        n - argmax_eq,
        n
    );
}
