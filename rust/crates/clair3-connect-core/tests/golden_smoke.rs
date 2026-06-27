use std::env;
use std::fs::File;
use std::io::{BufReader, Read, Seek};
use std::path::{Path, PathBuf};

use ndarray::Array2;
use ndarray_npy::ReadNpyExt;
use serde_json::Value;
use sha2::{Digest, Sha256};

fn sha256_file(path: &Path) -> std::io::Result<String> {
    let mut hasher = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    let mut fh = BufReader::new(File::open(path)?);
    loop {
        let n = fh.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(hex::encode(hasher.finalize()))
}

fn read_npy_header(path: &Path) -> std::io::Result<NpyHeader> {
    let mut fh = File::open(path)?;
    let mut magic = [0u8; 8];
    fh.read_exact(&mut magic)?;
    assert_eq!(
        &magic[..6],
        b"\x93NUMPY",
        "{}: not an .npy file",
        path.display()
    );
    let major = magic[6];
    let _minor = magic[7];

    let header_len = match major {
        1 => {
            let mut buf = [0u8; 2];
            fh.read_exact(&mut buf)?;
            u16::from_le_bytes(buf) as u64
        }
        2 | 3 => {
            let mut buf = [0u8; 4];
            fh.read_exact(&mut buf)?;
            u32::from_le_bytes(buf) as u64
        }
        other => panic!("{}: unknown .npy major version {}", path.display(), other),
    };

    let mut header = vec![0u8; header_len as usize];
    fh.read_exact(&mut header)?;
    let header_str = std::str::from_utf8(&header).expect("header utf-8");

    let descr = parse_str_field(header_str, "descr");
    let fortran_order = parse_bool_field(header_str, "fortran_order");
    let shape = parse_shape_field(header_str);
    let data_offset = fh.stream_position()?;

    Ok(NpyHeader {
        descr,
        fortran_order,
        shape,
        data_offset,
        major_version: major,
    })
}

#[derive(Debug)]
struct NpyHeader {
    descr: String,
    fortran_order: bool,
    shape: Vec<usize>,
    data_offset: u64,
    major_version: u8,
}

fn parse_str_field(header: &str, key: &str) -> String {

    let needle = format!("'{}':", key);
    let i = header
        .find(&needle)
        .unwrap_or_else(|| panic!("missing {} in {}", key, header));
    let rest = &header[i + needle.len()..];
    let q1 = rest.find('\'').unwrap();
    let q2 = rest[q1 + 1..].find('\'').unwrap() + q1 + 1;
    rest[q1 + 1..q2].to_string()
}

fn parse_bool_field(header: &str, key: &str) -> bool {
    let needle = format!("'{}':", key);
    let i = header
        .find(&needle)
        .unwrap_or_else(|| panic!("missing {} in {}", key, header));
    let rest = &header[i + needle.len()..].trim_start();
    rest.starts_with("True")
}

fn parse_shape_field(header: &str) -> Vec<usize> {
    let needle = "'shape':";
    let i = header
        .find(needle)
        .unwrap_or_else(|| panic!("missing shape in {}", header));
    let rest = &header[i + needle.len()..];
    let p1 = rest.find('(').unwrap();
    let p2 = rest.find(')').unwrap();
    let inner = &rest[p1 + 1..p2];
    inner
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<usize>().expect("shape int"))
        .collect()
}

struct FileEntry {
    path: String,
    sha256_file: String,
    sha256_bytes: Option<String>,
    nbytes: Option<u64>,

    dtype_shape: Option<(String, Vec<i64>)>,
}

fn sha256_hex(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    hex::encode(hasher.finalize())
}

fn collect_file_entries(v: &Value, out: &mut Vec<FileEntry>) {
    match v {
        Value::Object(map) => {
            if let (Some(Value::String(path)), Some(Value::String(sha))) =
                (map.get("path"), map.get("sha256_file"))
            {
                let dtype_shape = match (map.get("dtype"), map.get("shape")) {
                    (Some(Value::String(d)), Some(Value::Array(s))) => Some((
                        d.clone(),
                        s.iter().filter_map(|x| x.as_i64()).collect::<Vec<_>>(),
                    )),
                    _ => None,
                };
                let sha256_bytes = map
                    .get("sha256_bytes")
                    .and_then(|v| v.as_str())
                    .map(|s| s.to_string());
                let nbytes = map.get("nbytes").and_then(|v| v.as_u64());
                out.push(FileEntry {
                    path: path.clone(),
                    sha256_file: sha.clone(),
                    sha256_bytes,
                    nbytes,
                    dtype_shape,
                });
            }
            for child in map.values() {
                collect_file_entries(child, out);
            }
        }
        Value::Array(arr) => {
            for el in arr {
                collect_file_entries(el, out);
            }
        }
        _ => {}
    }
}

fn golden_dir() -> Option<PathBuf> {
    env::var_os("CLAIR3_CONNECT_GOLDEN_DIR").map(PathBuf::from)
}

#[test]
fn golden_artifacts_match_manifest() {
    let Some(dir) = golden_dir() else {
        eprintln!(
            "[skip] CLAIR3_CONNECT_GOLDEN_DIR not set. Run \
             `python tools/golden/produce_golden.py ... --output <dir>` then \
             `CLAIR3_CONNECT_GOLDEN_DIR=<dir> cargo test --workspace`."
        );
        return;
    };
    let manifest_path = dir.join("MANIFEST.json");
    let manifest_str = std::fs::read_to_string(&manifest_path)
        .unwrap_or_else(|e| panic!("read {}: {}", manifest_path.display(), e));
    let manifest: Value = serde_json::from_str(&manifest_str).expect("manifest parses");

    let schema = manifest
        .get("schema_version")
        .and_then(|v| v.as_str())
        .unwrap_or("");
    assert_eq!(schema, "phase0/0.1", "manifest schema version mismatch");

    let mut entries = Vec::new();
    collect_file_entries(&manifest, &mut entries);
    assert!(!entries.is_empty(), "manifest contained no file entries");

    let mut npy_entries_checked = 0usize;
    for entry in &entries {
        let abs = dir.join(&entry.path);
        let actual_sha = sha256_file(&abs).unwrap_or_else(|e| {
            panic!("sha256 of {} failed: {}", abs.display(), e);
        });
        assert_eq!(
            actual_sha,
            entry.sha256_file,
            "sha256 mismatch for {}: expected {}, got {}",
            abs.display(),
            entry.sha256_file,
            actual_sha
        );

        if abs.extension().map(|s| s == "npy").unwrap_or(false) {
            if let Some((manifest_dtype, manifest_shape)) = &entry.dtype_shape {
                let h = read_npy_header(&abs).unwrap();
                let header_shape: Vec<i64> = h.shape.iter().map(|&x| x as i64).collect();
                assert_eq!(
                    header_shape,
                    *manifest_shape,
                    "{}: .npy header shape {:?} != manifest shape {:?}",
                    abs.display(),
                    header_shape,
                    manifest_shape
                );
                let descr = &h.descr;
                let (item_bytes, ok_dtype) = match manifest_dtype.as_str() {
                    "int32" => (4u64, descr == "<i4" || descr == "=i4"),
                    "int64" => (8u64, descr == "<i8" || descr == "=i8"),
                    "float32" => (4u64, descr == "<f4" || descr == "=f4"),
                    "float64" => (8u64, descr == "<f8" || descr == "=f8"),
                    "int8" => (1u64, descr == "|i1" || descr == "<i1"),
                    other => {
                        eprintln!(
                            "[warn] {}: unhandled manifest dtype {}",
                            abs.display(),
                            other
                        );
                        (0u64, true)
                    }
                };
                assert!(
                    ok_dtype,
                    "{}: .npy descr {} not compatible with manifest dtype {}",
                    abs.display(),
                    descr,
                    manifest_dtype
                );
                assert!(
                    !h.fortran_order,
                    "{}: unexpected fortran-order .npy",
                    abs.display()
                );
                assert!(
                    h.data_offset >= 10,
                    "{}: data offset too small",
                    abs.display()
                );
                let _ = h.major_version;

                if item_bytes > 0 {
                    let n_elements: u64 = manifest_shape.iter().map(|&d| d as u64).product();
                    let expected_payload_bytes = n_elements * item_bytes;
                    if let Some(manifest_nbytes) = entry.nbytes {
                        assert_eq!(
                            manifest_nbytes,
                            expected_payload_bytes,
                            "{}: manifest nbytes {} != shape*itemsize {}",
                            abs.display(),
                            manifest_nbytes,
                            expected_payload_bytes
                        );
                    }

                    let mut payload = Vec::with_capacity(expected_payload_bytes as usize);
                    let mut fh = File::open(&abs).unwrap();
                    fh.seek(std::io::SeekFrom::Start(h.data_offset)).unwrap();
                    fh.read_to_end(&mut payload).unwrap();
                    assert_eq!(
                        payload.len() as u64,
                        expected_payload_bytes,
                        "{}: actual payload {} != expected {}",
                        abs.display(),
                        payload.len(),
                        expected_payload_bytes
                    );

                    if let Some(manifest_sha_bytes) = &entry.sha256_bytes {
                        let actual = sha256_hex(&payload);
                        assert_eq!(
                            &actual,
                            manifest_sha_bytes,
                            "{}: sha256_bytes mismatch (got {} expected {})",
                            abs.display(),
                            actual,
                            manifest_sha_bytes
                        );
                    }
                }
                npy_entries_checked += 1;
            }
        }
    }

    eprintln!(
        "[ok] verified {} files ({} .npy headers) in {}",
        entries.len(),
        npy_entries_checked,
        dir.display()
    );
    assert!(
        npy_entries_checked >= 4,
        "expected at least 4 .npy header checks (raw matrix, major, minor, X)"
    );
}

#[test]
fn pileup_raw_matrix_has_negative_cells() {
    let Some(dir) = golden_dir() else {
        eprintln!("[skip] CLAIR3_CONNECT_GOLDEN_DIR not set");
        return;
    };
    let path = dir.join("pileup_raw_chunk_0000_matrix.npy");
    if !path.exists() {
        eprintln!("[skip] {} not present", path.display());
        return;
    }
    let m = ndarray::Array2::<i64>::read_npy(File::open(&path).unwrap())
        .expect("read raw matrix as i64");
    let total = m.len() as f64;
    let neg = m.iter().filter(|&&v| v < 0).count() as f64;
    let frac = neg / total;
    eprintln!(
        "[info] raw matrix: shape={:?} total_cells={} neg_cells={} (frac={:.4})",
        m.shape(),
        total as usize,
        neg as usize,
        frac
    );
    assert!(
        frac > 0.05,
        "expected at least 5% negative cells from the -1*forward_sum wrap; got {:.4}",
        frac
    );

    let neg_cols: Vec<usize> = (0..m.ncols())
        .filter(|&c| m.column(c).iter().any(|&v| v < 0))
        .collect();
    eprintln!("[info] columns with negative cells: {:?}", neg_cols);
    for &c in &neg_cols {
        assert!(
            c < 4 || (9..13).contains(&c),
            "negative cells in unexpected column {} (should be 0..4 fwd ACGT or 9..13 rev acgt)",
            c
        );
    }
}

#[test]
fn pileup_y_is_valid_softmax_output() {
    let Some(dir) = golden_dir() else {
        eprintln!("[skip] CLAIR3_CONNECT_GOLDEN_DIR not set");
        return;
    };
    let manifest: Value =
        serde_json::from_str(&std::fs::read_to_string(dir.join("MANIFEST.json")).unwrap()).unwrap();

    let y_meta = manifest["stages"]["pileup_predict"]["data"]["Y_cpu"].clone();
    let y_path = dir.join(y_meta["path"].as_str().expect("Y_cpu path"));
    let n_candidates = y_meta["shape"][0].as_i64().unwrap() as usize;
    let n_cols = y_meta["shape"][1].as_i64().unwrap() as usize;
    assert_eq!(
        n_cols, 24,
        "pileup output should be GT21(21) + GENOTYPE(3) = 24"
    );

    let y = Array2::<f32>::read_npy(File::open(&y_path).unwrap())
        .expect("read Y_pileup_cpu.npy as f32");
    assert_eq!(y.shape(), &[n_candidates, n_cols], "Y shape mismatch");

    let eps = 1e-3_f32;
    let mut bad_rows: Vec<(usize, f32, f32)> = Vec::new();
    for i in 0..n_candidates {
        let gt21_sum: f32 = y.row(i).iter().take(21).sum();
        let gt_sum: f32 = y.row(i).iter().skip(21).take(3).sum();
        if (gt21_sum - 1.0).abs() > eps || (gt_sum - 1.0).abs() > eps {
            bad_rows.push((i, gt21_sum, gt_sum));
        }
    }
    if !bad_rows.is_empty() {
        for &(i, g21, g) in bad_rows.iter().take(5) {
            eprintln!("non-softmax row {}: gt21_sum={}, gt_sum={}", i, g21, g);
        }
    }
    assert!(
        bad_rows.is_empty(),
        "{}/{} rows have non-softmax head sums (tol {:.0e})",
        bad_rows.len(),
        n_candidates,
        eps
    );
    eprintln!(
        "[ok] Y_pileup_cpu.npy: shape {:?}, {} rows softmax-valid",
        y.shape(),
        n_candidates
    );
}
