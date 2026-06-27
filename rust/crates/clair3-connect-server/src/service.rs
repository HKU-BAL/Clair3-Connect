use crate::backpressure::{BackpressureConfig, Gate};
use arc_swap::ArcSwap;
use clair3_connect_infer::batcher::{BatchConfig, BatchError, BatchedWorker};
use clair3_connect_infer::model_registry::{ModelEntry, Registry};
use clair3_connect_infer::torchscript::LoadedModel;
use clair3_connect_infer::worker::Worker;
use clair3_connect_transport::aead::{self, AeadKey};
use clair3_connect_transport::checksum::sha256_hex;
use clair3_connect_transport::proto::{
    clair3_connect_inference_server::Clair3ConnectInference, HealthRequest, HealthResponse, PredictRequest,
    PredictResponse, PredictStatus, TensorDType,
};
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;
use tch::{Device, Kind, Tensor};
use tonic::{Request, Response, Status};

fn device_label(d: Device) -> String {
    match d {
        Device::Cpu => "cpu".to_string(),
        Device::Cuda(n) => format!("cuda:{n}"),

        other => format!("{other:?}").to_lowercase(),
    }
}

struct LoadedEntry {

    _model_id: String,

    batched_per_device: Vec<BatchedWorker>,

    devices: Vec<Device>,

    per_device_gates: Vec<Gate>,

    rr_counter: std::sync::Arc<std::sync::atomic::AtomicUsize>,

    model_sha256: String,
}

impl LoadedEntry {

    #[allow(clippy::result_large_err)]
    fn pick(
        &self,
    ) -> Result<
        (
            usize,
            &BatchedWorker,
            Device,
            tokio::sync::OwnedSemaphorePermit,
        ),
        Status,
    > {
        let n = self.batched_per_device.len();
        let start = if n == 1 {
            0
        } else {
            self.rr_counter
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed)
                % n
        };
        for offset in 0..n {
            let idx = (start + offset) % n;
            if let Ok(permit) = self.per_device_gates[idx].try_acquire() {
                return Ok((
                    idx,
                    &self.batched_per_device[idx],
                    self.devices[idx],
                    permit,
                ));
            }
        }
        Err(Status::resource_exhausted(format!(
            "all {n} per-device gates at capacity"
        )))
    }

}

#[derive(Clone)]
pub struct PredictService {

    entries: Arc<ArcSwap<HashMap<String, Arc<LoadedEntry>>>>,
    gate: Gate,

    aead_key: Option<Arc<AeadKey>>,

    devices: Arc<Vec<Device>>,
    bp: BackpressureConfig,
    batch: BatchConfig,
}

#[derive(Debug, Default)]
pub struct ReloadSummary {

    pub added: Vec<String>,

    pub updated: Vec<String>,

    pub removed: Vec<String>,

    pub unchanged: usize,
}

fn load_one_entry(
    registry: &Registry,
    entry: &ModelEntry,
    registry_yaml_dir: &Path,
    devices: &[Device],
    bp: &BackpressureConfig,
    batch: &BatchConfig,
) -> anyhow::Result<Option<Arc<LoadedEntry>>> {
    let candidate_ts = registry_yaml_dir.join(&entry.torchscript.path);
    if !candidate_ts.exists() {
        tracing::warn!(
            model_id = %entry.model_id,
            path = %candidate_ts.display(),
            "skipping registry entry — .ts file does not exist (template entry?)"
        );
        return Ok(None);
    }
    let ts_path = registry
        .verify_torchscript(entry, registry_yaml_dir)
        .map_err(|e| anyhow::anyhow!("verify {}: {}", entry.model_id, e))?;
    let model_sha256 = entry.torchscript.sha256.clone();
    let mut batched_per_device = Vec::with_capacity(devices.len());
    let mut per_device_gates = Vec::with_capacity(devices.len());
    for &device in devices {
        let label = device_label(device);
        tracing::info!(
            model_id = %entry.model_id,
            path = %ts_path.display(),
            device = %label,
            "loading TorchScript"
        );
        let model = LoadedModel::load(&ts_path, entry.clone(), device)
            .map_err(|e| anyhow::anyhow!("LoadedModel::load {} on {}: {}", entry.model_id, label, e))?;
        let worker = Worker::new(model);
        batched_per_device.push(BatchedWorker::spawn(worker, batch.clone()));
        per_device_gates.push(Gate::new(bp.clone()));
    }
    Ok(Some(Arc::new(LoadedEntry {
        _model_id: entry.model_id.clone(),
        batched_per_device,
        devices: devices.to_vec(),
        per_device_gates,
        rr_counter: std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0)),
        model_sha256,
    })))
}

impl PredictService {

    pub fn load(
        registry: &Registry,
        registry_yaml_dir: &Path,
        device: Device,
    ) -> anyhow::Result<Self> {
        Self::load_with_backpressure(
            registry,
            registry_yaml_dir,
            device,
            BackpressureConfig::default(),
        )
    }

    pub fn load_with_backpressure(
        registry: &Registry,
        registry_yaml_dir: &Path,
        device: Device,
        bp: BackpressureConfig,
    ) -> anyhow::Result<Self> {
        Self::load_with_options(
            registry,
            registry_yaml_dir,
            &[device],
            bp,
            BatchConfig::default(),
        )
    }

    pub fn load_with_options(
        registry: &Registry,
        registry_yaml_dir: &Path,
        devices: &[Device],
        bp: BackpressureConfig,
        batch: BatchConfig,
    ) -> anyhow::Result<Self> {
        if devices.is_empty() {
            anyhow::bail!("PredictService::load_with_options requires at least one device");
        }

        let cuda_count = clair3_connect_infer::_ensure_cuda_linked();
        tracing::debug!(cuda_count, "CUDA link probe");
        let mut map: HashMap<String, Arc<LoadedEntry>> = HashMap::new();
        for entry in registry.models.iter() {

            if let Some(loaded) = load_one_entry(registry, entry, registry_yaml_dir, devices, &bp, &batch)? {
                map.insert(entry.model_id.clone(), loaded);
            }
        }
        Ok(Self {
            entries: Arc::new(ArcSwap::from_pointee(map)),
            gate: Gate::new(bp.clone()),
            aead_key: None,
            devices: Arc::new(devices.to_vec()),
            bp,
            batch,
        })
    }

    pub fn reload_from_registry_file(&self, registry_path: &Path) -> anyhow::Result<ReloadSummary> {
        let registry = Registry::load(registry_path)
            .map_err(|e| anyhow::anyhow!("load registry {}: {}", registry_path.display(), e))?;
        let dir = registry_path.parent().unwrap_or_else(|| Path::new("."));
        let current = self.entries.load();
        let mut new_map: HashMap<String, Arc<LoadedEntry>> = HashMap::new();
        let mut summary = ReloadSummary::default();
        for entry in registry.models.iter() {

            if let Some(existing) = current.get(&entry.model_id) {
                if existing.model_sha256 == entry.torchscript.sha256 {
                    new_map.insert(entry.model_id.clone(), existing.clone());
                    summary.unchanged += 1;
                    continue;
                }
            }
            let was_present = current.contains_key(&entry.model_id);
            if let Some(loaded) =
                load_one_entry(&registry, entry, dir, &self.devices, &self.bp, &self.batch)?
            {
                new_map.insert(entry.model_id.clone(), loaded);
                if was_present {
                    summary.updated.push(entry.model_id.clone());
                } else {
                    summary.added.push(entry.model_id.clone());
                }
            }
        }
        summary.removed = current
            .keys()
            .filter(|k| !new_map.contains_key(*k))
            .cloned()
            .collect();
        self.entries.store(Arc::new(new_map));
        Ok(summary)
    }

    pub fn set_aead_key(&mut self, key: AeadKey) {
        self.aead_key = Some(Arc::new(key));
    }

    pub fn aead_enabled(&self) -> bool {
        self.aead_key.is_some()
    }

    pub fn max_inflight(&self) -> usize {
        self.gate.max_inflight()
    }

    pub fn metrics_snapshot(&self) -> crate::metrics_server::MetricsSnapshot {
        use std::sync::atomic::Ordering;
        let entries = self.entries.load();
        let mut per_model = Vec::with_capacity(entries.len());
        for (id, entry) in entries.iter() {
            for (i, worker) in entry.batched_per_device.iter().enumerate() {
                let m = worker.metrics();
                let gate = &entry.per_device_gates[i];
                per_model.push(crate::metrics_server::ModelMetrics {
                    model_id: id.clone(),
                    model_sha256: entry.model_sha256.clone(),
                    device: device_label(entry.devices[i]),
                    batches: m.batches.load(Ordering::Relaxed),
                    items: m.items.load(Ordering::Relaxed),
                    compute_ns: m.compute_ns.load(Ordering::Relaxed),
                    gate_max_inflight: gate.max_inflight(),
                    gate_available_permits: gate.available_permits(),
                });
            }
        }

        per_model.sort_by(|a, b| {
            a.model_id
                .cmp(&b.model_id)
                .then_with(|| a.device.cmp(&b.device))
        });
        crate::metrics_server::MetricsSnapshot {
            gate_max_inflight: self.gate.max_inflight(),
            gate_available_permits: self.gate.available_permits(),
            per_model,
        }
    }
}

fn dtype_to_kind(dtype: i32) -> Option<Kind> {
    let d = TensorDType::try_from(dtype).ok()?;
    match d {
        TensorDType::TensorDtypeInt8 => Some(Kind::Int8),
        TensorDType::TensorDtypeInt16 => Some(Kind::Int16),
        TensorDType::TensorDtypeInt32 => Some(Kind::Int),
        TensorDType::TensorDtypeFloat32 => Some(Kind::Float),
        TensorDType::TensorDtypeUnspecified => None,
    }
}

fn kind_to_dtype(kind: Kind) -> i32 {
    match kind {
        Kind::Int8 => TensorDType::TensorDtypeInt8 as i32,
        Kind::Int16 => TensorDType::TensorDtypeInt16 as i32,
        Kind::Int => TensorDType::TensorDtypeInt32 as i32,
        Kind::Float => TensorDType::TensorDtypeFloat32 as i32,
        _ => TensorDType::TensorDtypeUnspecified as i32,
    }
}

fn decode_tensor(dtype: i32, shape: &[i64], payload: &[u8], device: Device) -> Result<Tensor, String> {
    let kind = dtype_to_kind(dtype).ok_or_else(|| "unspecified dtype".to_string())?;
    let elt_size = match kind {
        Kind::Int8 => 1,
        Kind::Int16 => 2,
        Kind::Int => 4,
        Kind::Float => 4,
        _ => return Err("unsupported dtype".to_string()),
    };
    let n_elements: usize = {
        let mut n: usize = 1;
        for d in shape {
            if *d < 0 {
                return Err("negative dim".to_string());
            }
            n = n
                .checked_mul(*d as usize)
                .ok_or_else(|| "shape overflow".to_string())?;
        }
        n
    };
    let expected_bytes = n_elements
        .checked_mul(elt_size)
        .ok_or_else(|| "shape * elt_size overflow".to_string())?;
    if payload.len() != expected_bytes {
        return Err(format!(
            "payload {} bytes, expected {}",
            payload.len(),
            expected_bytes
        ));
    }

    let tensor = match kind {
        Kind::Int8 => {
            let v: Vec<i8> = payload.iter().map(|b| *b as i8).collect();
            Tensor::from_slice(&v).reshape(shape).to_device(device)
        }
        Kind::Int16 => {
            let v: Vec<i16> = payload
                .chunks_exact(2)
                .map(|c| i16::from_le_bytes([c[0], c[1]]))
                .collect();
            Tensor::from_slice(&v).reshape(shape).to_device(device)
        }
        Kind::Int => {
            let v: Vec<i32> = payload
                .chunks_exact(4)
                .map(|c| i32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                .collect();
            Tensor::from_slice(&v).reshape(shape).to_device(device)
        }
        Kind::Float => {
            let v: Vec<f32> = payload
                .chunks_exact(4)
                .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                .collect();
            Tensor::from_slice(&v).reshape(shape).to_device(device)
        }
        _ => unreachable!(),
    };
    Ok(tensor)
}

fn encode_tensor(y: &Tensor) -> Result<(i32, Vec<i64>, Vec<u8>), String> {
    let kind = y.kind();
    let shape = y.size();
    let cpu = y.to_kind(kind).to_device(Device::Cpu).contiguous();
    let payload = match kind {
        Kind::Float => {
            let n: i64 = shape.iter().copied().product();
            let v: Vec<f32> = Vec::<f32>::try_from(cpu.reshape([n]))
                .map_err(|e| format!("tch to Vec<f32>: {e}"))?;
            let mut out = Vec::with_capacity(v.len() * 4);
            for f in &v {
                out.extend_from_slice(&f.to_le_bytes());
            }
            out
        }
        Kind::Int8 => {
            let n: i64 = shape.iter().copied().product();
            let v: Vec<i8> = Vec::<i8>::try_from(cpu.reshape([n]))
                .map_err(|e| format!("tch to Vec<i8>: {e}"))?;
            v.iter().map(|i| *i as u8).collect()
        }
        Kind::Int => {
            let n: i64 = shape.iter().copied().product();
            let v: Vec<i32> = Vec::<i32>::try_from(cpu.reshape([n]))
                .map_err(|e| format!("tch to Vec<i32>: {e}"))?;
            let mut out = Vec::with_capacity(v.len() * 4);
            for i in &v {
                out.extend_from_slice(&i.to_le_bytes());
            }
            out
        }
        Kind::Int16 => {
            let n: i64 = shape.iter().copied().product();
            let v: Vec<i16> = Vec::<i16>::try_from(cpu.reshape([n]))
                .map_err(|e| format!("tch to Vec<i16>: {e}"))?;
            let mut out = Vec::with_capacity(v.len() * 2);
            for i in &v {
                out.extend_from_slice(&i.to_le_bytes());
            }
            out
        }
        _ => return Err(format!("unsupported output kind {:?}", kind)),
    };
    Ok((kind_to_dtype(kind), shape, payload))
}

#[tonic::async_trait]
impl Clair3ConnectInference for PredictService {
    async fn predict(
        &self,
        request: Request<PredictRequest>,
    ) -> Result<Response<PredictResponse>, Status> {

        let _permit = self.gate.try_acquire()?;

        let mut req = request.into_inner();

        let entry = self
            .entries
            .load()
            .get(&req.model_id)
            .cloned()
            .ok_or_else(|| Status::not_found(format!("unknown model_id {:?}", req.model_id)))?;

        let client_encrypted = !req.aead_nonce.is_empty();
        match (self.aead_key.as_ref(), client_encrypted) {
            (Some(key), true) => {
                let aad = aead::request_aad(
                    &req.model_id,
                    req.dtype,
                    &req.shape,
                    &req.payload_sha256,
                );
                let plaintext = match aead::open(key, &req.payload, &req.aead_nonce, &aad)
                {
                    Ok(p) => p,
                    Err(e) => {
                        return Ok(Response::new(reject_response(
                            &req.request_id,
                            "aead_open_failed",
                            &format!("{e}"),
                            "",
                            &entry.model_sha256,
                        )));
                    }
                };
                req.payload = plaintext.into();
                req.aead_nonce = Default::default();
            }
            (Some(_), false) => {
                return Ok(Response::new(reject_response(
                    &req.request_id,
                    "aead_required",
                    "server requires AEAD; client sent plaintext",
                    "",
                    &entry.model_sha256,
                )));
            }
            (None, true) => {
                return Ok(Response::new(reject_response(
                    &req.request_id,
                    "aead_not_supported",
                    "client sent AEAD ciphertext; server has no key",
                    "",
                    &entry.model_sha256,
                )));
            }
            (None, false) => {  }
        }

        let actual_sha = sha256_hex(&req.payload);
        if !req.payload_sha256.is_empty() && req.payload_sha256 != actual_sha {
            return Ok(Response::new(reject_response(
                &req.request_id,
                "payload_sha_mismatch",
                &format!("expected {} got {}", req.payload_sha256, actual_sha),
                &actual_sha,
                &entry.model_sha256,
            )));
        }

        let (_dev_idx, worker, device, _dev_permit) = match entry.pick() {
            Ok(t) => t,
            Err(status) => return Err(status),
        };
        let x = match decode_tensor(req.dtype, &req.shape, &req.payload, device) {
            Ok(t) => t,
            Err(e) => {
                return Ok(Response::new(reject_response(
                    &req.request_id,
                    "decode_failed",
                    &e,
                    &actual_sha,
                    &entry.model_sha256,
                )));
            }
        };

        let y = match worker.predict(x).await {
            Ok(t) => t,
            Err(BatchError::Full) => {
                return Err(Status::resource_exhausted(
                    "batcher queue full — retry shortly",
                ));
            }
            Err(BatchError::ShutDown) => {
                return Ok(Response::new(failed_response(
                    &req.request_id,
                    "forward_failed",
                    "batcher shut down",
                    &actual_sha,
                    &entry.model_sha256,
                )));
            }
            Err(BatchError::Forward(msg)) => {
                return Ok(Response::new(failed_response(
                    &req.request_id,
                    "forward_failed",
                    &msg,
                    &actual_sha,
                    &entry.model_sha256,
                )));
            }
        };

        let (y_dtype, y_shape, y_payload) = match encode_tensor(&y) {
            Ok(t) => t,
            Err(e) => {
                return Ok(Response::new(failed_response(
                    &req.request_id,
                    "encode_failed",
                    &e,
                    &actual_sha,
                    &entry.model_sha256,
                )));
            }
        };
        let y_sha = sha256_hex(&y_payload);

        let (y_payload_wire, y_nonce_wire): (Vec<u8>, Vec<u8>) = match self.aead_key.as_ref() {
            Some(key) => {
                let aad = aead::response_aad(&req.request_id, y_dtype, &y_shape, &y_sha);
                let (ct, nonce) = aead::seal(key, &y_payload, &aad);
                (ct, nonce)
            }
            None => (y_payload, Vec::new()),
        };

        Ok(Response::new(PredictResponse {
            request_id: req.request_id,
            status: PredictStatus::Ok as i32,
            error_code: String::new(),
            error_message: String::new(),
            y_dtype,
            y_shape,
            y_payload: y_payload_wire.into(),
            y_sha256: y_sha,
            input_payload_sha256: actual_sha,
            model_sha256: entry.model_sha256.clone(),
            y_aead_nonce: y_nonce_wire.into(),
        }))
    }

    async fn health(
        &self,
        _request: Request<HealthRequest>,
    ) -> Result<Response<HealthResponse>, Status> {
        let mut registered_models = HashMap::new();
        for (id, e) in self.entries.load().iter() {
            registered_models.insert(id.clone(), e.model_sha256.clone());
        }
        Ok(Response::new(HealthResponse {
            server_version: env!("CARGO_PKG_VERSION").to_string(),
            registered_models,
        }))
    }
}

fn reject_response(
    req_id: &str,
    code: &str,
    msg: &str,
    input_sha: &str,
    model_sha: &str,
) -> PredictResponse {
    PredictResponse {
        request_id: req_id.to_string(),
        status: PredictStatus::Rejected as i32,
        error_code: code.to_string(),
        error_message: msg.to_string(),
        y_dtype: TensorDType::TensorDtypeUnspecified as i32,
        y_shape: Vec::new(),
        y_payload: Vec::new().into(),
        y_sha256: String::new(),
        input_payload_sha256: input_sha.to_string(),
        model_sha256: model_sha.to_string(),
        y_aead_nonce: Default::default(),
    }
}

fn failed_response(
    req_id: &str,
    code: &str,
    msg: &str,
    input_sha: &str,
    model_sha: &str,
) -> PredictResponse {
    PredictResponse {
        request_id: req_id.to_string(),
        status: PredictStatus::Failed as i32,
        error_code: code.to_string(),
        error_message: msg.to_string(),
        y_dtype: TensorDType::TensorDtypeUnspecified as i32,
        y_shape: Vec::new(),
        y_payload: Vec::new().into(),
        y_sha256: String::new(),
        input_payload_sha256: input_sha.to_string(),
        model_sha256: model_sha.to_string(),
        y_aead_nonce: Default::default(),
    }
}
