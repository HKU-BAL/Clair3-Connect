use std::path::PathBuf;
use tch::{Device, Kind, Tensor};

use clair3_connect_infer::model_registry::Registry;
use clair3_connect_infer::torchscript::LoadedModel;
use clair3_connect_infer::worker::Worker;

const MODEL_ID: &str = "full_alignment_ont_hac_v500";

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

fn artifacts() -> PathBuf {
    workspace_root().join("test_data/ont_chr20_demo/artifacts")
}

fn registry_yaml() -> PathBuf {
    workspace_root().join("models/model_registry.example.yaml")
}

fn read_npy_f32(path: &std::path::Path, expected_cols: usize) -> (usize, Vec<f32>) {
    let bytes = std::fs::read(path).expect("read npy");
    let major = bytes[6];
    let header_len = if major == 1 {
        u16::from_le_bytes([bytes[8], bytes[9]]) as usize + 10
    } else {
        u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize + 12
    };
    let data = &bytes[header_len..];
    let total = data.len() / 4;
    let rows = total / expected_cols;
    let mut out = Vec::with_capacity(total);
    for i in 0..total {
        let off = i * 4;
        out.push(f32::from_le_bytes([
            data[off],
            data[off + 1],
            data[off + 2],
            data[off + 3],
        ]));
    }
    (rows, out)
}

fn read_npy_i8(path: &std::path::Path) -> Vec<i8> {
    let bytes = std::fs::read(path).expect("read npy");
    let major = bytes[6];
    let header_len = if major == 1 {
        u16::from_le_bytes([bytes[8], bytes[9]]) as usize + 10
    } else {
        u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize + 12
    };
    bytes[header_len..].iter().map(|b| *b as i8).collect()
}

#[test]
fn g7_fa_y_parity_cpu() {
    let arts = artifacts();
    let yaml = registry_yaml();
    let matrix_path = arts.join("fa_haplotag_matrix.npy");
    let y_path = arts.join("Y_full_alignment_cpu.npy");
    if !matrix_path.exists() || !y_path.exists() || !yaml.exists() {
        eprintln!("Skipping — fixture incomplete");
        return;
    }

    let reg = Registry::load(&yaml).expect("registry");
    let entry = reg.get(MODEL_ID).expect("FA model registered");
    let ts = match reg.verify_torchscript(entry, yaml.parent().unwrap()) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("Skipping — TorchScript verify failed: {e}");
            return;
        }
    };
    let model = LoadedModel::load(&ts, entry.clone(), Device::Cpu).expect("CModule load");
    let worker = Worker::new(model);

    let x_bytes = read_npy_i8(&matrix_path);
    let n = 64;
    let d = 89;
    let p = 33;
    let c = 8;
    assert_eq!(x_bytes.len(), n * d * p * c, "matrix size unexpected");
    let x = Tensor::from_slice(&x_bytes)
        .reshape([n as i64, d as i64, p as i64, c as i64])
        .to_kind(Kind::Int8);
    let y_rust = worker.forward_batch(&x).expect("FA forward");
    assert_eq!(y_rust.size(), vec![n as i64, 90]);

    let y_rust_flat = y_rust
        .to_kind(Kind::Float)
        .contiguous()
        .reshape([n as i64 * 90]);
    let y_rust_v: Vec<f32> = Vec::<f32>::try_from(&y_rust_flat).expect("rust Y to Vec");

    let (rows, y_py) = read_npy_f32(&y_path, 90);
    assert_eq!(rows, n, "python Y rows differ");
    assert_eq!(y_rust_v.len(), y_py.len());

    let mut max_abs = 0.0_f32;
    let mut arg_match = 0;
    for r in 0..n {
        let r_off = r * 90;
        let mut arg_rust = 0;
        let mut arg_py = 0;
        for (i, (a, b)) in y_rust_v[r_off..r_off + 90]
            .iter()
            .zip(y_py[r_off..r_off + 90].iter())
            .enumerate()
        {
            assert!(a.is_finite() && b.is_finite(), "non-finite at row {r} col {i}");
            let d = (a - b).abs();
            if d > max_abs {
                max_abs = d;
            }
            if y_rust_v[r_off + i] > y_rust_v[r_off + arg_rust] {
                arg_rust = i;
            }
            if y_py[r_off + i] > y_py[r_off + arg_py] {
                arg_py = i;
            }
        }
        if arg_rust == arg_py {
            arg_match += 1;
        }
    }
    eprintln!("G7 max_abs_diff = {max_abs:e}; argmax match = {arg_match}/{n}");
    assert!(
        max_abs < 1e-5,
        "G7 FA Y parity exceeded CPU tolerance: max_abs={max_abs:e}"
    );
}
