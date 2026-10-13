use crate::checksum::sha256_hex;
use crate::codec;
use crate::proto::clair3_connect_inference_client::Clair3ConnectInferenceClient;
use crate::proto::{HealthRequest, PredictRequest, PredictResponse, PredictStatus, TensorDType};
use std::time::Duration;
use thiserror::Error;
use tonic::transport::Channel;
use tonic::Code;

#[derive(Debug, Clone, Copy)]
pub struct RetryPolicy {

    pub max_retries: u32,

    pub initial_delay: Duration,

    pub max_delay: Duration,

    pub per_attempt_timeout: Option<Duration>,
}

impl Default for RetryPolicy {

    fn default() -> Self {
        Self {
            max_retries: 3,
            initial_delay: Duration::from_millis(100),
            max_delay: Duration::from_secs(5),
            per_attempt_timeout: Some(Duration::from_secs(30)),
        }
    }
}

impl RetryPolicy {

    pub const fn none() -> Self {
        Self {
            max_retries: 0,
            initial_delay: Duration::from_millis(0),
            max_delay: Duration::from_millis(0),
            per_attempt_timeout: None,
        }
    }

    pub fn next_backoff_for(&self, attempt: u32) -> Duration {
        let base = self.base_delay_for(attempt);
        let now_nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.subsec_nanos())
            .unwrap_or(0);
        let init = self.initial_delay.as_nanos() as u64;
        let jitter_nanos = if init == 0 { 0 } else { (now_nanos as u64) % init };
        let with_jitter = base + Duration::from_nanos(jitter_nanos);
        with_jitter.min(self.max_delay)
    }

    pub fn is_retryable_grpc(code: Code) -> bool {
        matches!(
            code,
            Code::Unavailable | Code::Aborted | Code::ResourceExhausted
        )
    }

    pub fn is_retryable_status(status: PredictStatus) -> bool {
        matches!(status, PredictStatus::Retry)
    }

    pub fn base_delay_for(self, attempt: u32) -> Duration {
        let scaled = self
            .initial_delay
            .saturating_mul(1_u32.checked_shl(attempt).unwrap_or(u32::MAX));
        scaled.min(self.max_delay)
    }
}

#[derive(Clone)]
pub struct PredictClient {
    inner: Clair3ConnectInferenceClient<Channel>,
    retry_policy: RetryPolicy,

    aead_key: Option<std::sync::Arc<crate::aead::AeadKey>>,
}

#[derive(Debug, Error)]
pub enum PredictError {

    #[error("transport: {0}")]
    Transport(#[from] tonic::transport::Error),

    #[error("grpc status: {0}")]
    Grpc(#[from] tonic::Status),

    #[error("response sha256 mismatch: declared {declared:?}, computed {computed:?}")]
    Sha256Mismatch {

        declared: String,

        computed: String,
    },

    #[error("server model_sha256 echo mismatch: expected {expected:?}, got {got:?}")]
    ModelShaMismatch {

        expected: String,

        got: String,
    },

    #[error("predict status {status:?}: {code} {message}")]
    NonOk {

        status: PredictStatus,

        code: String,

        message: String,
    },
}

fn configure_client_compression(
    client: Clair3ConnectInferenceClient<Channel>,
) -> Clair3ConnectInferenceClient<Channel> {
    use tonic::codec::CompressionEncoding;
    let client = client
        .accept_compressed(CompressionEncoding::Zstd)
        .accept_compressed(CompressionEncoding::Gzip);
    match std::env::var("CLAIR3_CONNECT_COMPRESS").ok().as_deref() {
        Some("none") | Some("off") | Some("0") => client,
        Some("gzip") => client.send_compressed(CompressionEncoding::Gzip),
        _ => client.send_compressed(CompressionEncoding::Zstd),
    }
}

impl PredictClient {

    pub async fn connect(endpoint: String) -> Result<Self, PredictError> {
        Self::connect_with_retry(endpoint, RetryPolicy::default()).await
    }

    pub async fn connect_with_retry(
        endpoint: String,
        retry_policy: RetryPolicy,
    ) -> Result<Self, PredictError> {
        let ep = tonic::transport::Endpoint::from_shared(endpoint).map_err(|e| {
            tonic::Status::invalid_argument(format!("bad endpoint: {e}"))
        })?

        .tcp_nodelay(true);
        let chan = connect_endpoint_with_retry(ep, retry_policy).await?;
        Ok(Self {
            inner: configure_client_compression(Clair3ConnectInferenceClient::new(chan)),
            retry_policy,
            aead_key: None,
        })
    }

    pub async fn connect_tls_with_retry(
        endpoint: String,
        tls: &crate::tls::ClientTlsOptions,
        retry_policy: RetryPolicy,
    ) -> Result<Self, PredictError> {
        let ep = tonic::transport::Endpoint::from_shared(endpoint)
            .map_err(|e| {
                tonic::Status::invalid_argument(format!("bad endpoint: {e}"))
            })?
            .tcp_nodelay(true)
            .tls_config(tls.build())
            .map_err(|e| {
                tonic::Status::invalid_argument(format!("tls_config: {e}"))
            })?;
        let chan = connect_endpoint_with_retry(ep, retry_policy).await?;
        Ok(Self {
            inner: configure_client_compression(Clair3ConnectInferenceClient::new(chan)),
            retry_policy,
            aead_key: None,
        })
    }

    pub fn from_channel(channel: Channel, retry_policy: RetryPolicy) -> Self {
        Self {
            inner: configure_client_compression(Clair3ConnectInferenceClient::new(channel)),
            retry_policy,
            aead_key: None,
        }
    }

    pub fn set_retry_policy(&mut self, policy: RetryPolicy) {
        self.retry_policy = policy;
    }

    pub fn set_aead_key(&mut self, key: crate::aead::AeadKey) {
        self.aead_key = Some(std::sync::Arc::new(key));
    }

    pub fn aead_enabled(&self) -> bool {
        self.aead_key.is_some()
    }

    pub async fn predict(
        &mut self,
        req: PredictRequest,
    ) -> Result<PredictResponse, PredictError> {
        self.predict_with_expected_model(req, None).await
    }

    pub async fn predict_with_expected_model(
        &mut self,
        mut req: PredictRequest,
        expected_model_sha256: Option<&str>,
    ) -> Result<PredictResponse, PredictError> {
        if req.payload_sha256.is_empty() {
            req.payload_sha256 = sha256_hex(&req.payload);
        }

        if let Some(key) = self.aead_key.as_ref() {
            let aad = crate::aead::request_aad(
                &req.model_id,
                req.dtype,
                &req.shape,
                &req.payload_sha256,
            );
            if crate::codec::request_zstd_enabled() {
                if let Ok(c) = crate::codec::compress_payload(&req.payload) {
                    if c.len() < req.payload.len() {
                        req.payload = c.into();
                        req.payload_encoding = crate::codec::PAYLOAD_ENCODING_ZSTD.to_string();
                    }
                }
            }
            let (ciphertext, nonce) = crate::aead::seal(key, &req.payload, &aad);
            req.payload = ciphertext.into();
            req.aead_nonce = nonce.into();
        }

        let mut last_err: Option<PredictError> = None;
        for attempt in 0..=self.retry_policy.max_retries {
            if attempt > 0 {
                let delay = self.next_backoff(attempt - 1);
                tokio::time::sleep(delay).await;
            }

            let timeout = self.retry_policy.per_attempt_timeout;
            let attempt_fut = self
                .try_predict_once(req.clone(), expected_model_sha256);
            let result = match timeout {
                Some(t) => match tokio::time::timeout(t, attempt_fut).await {
                    Ok(r) => r,
                    Err(_) => Err(PredictError::Grpc(tonic::Status::deadline_exceeded(
                        format!("predict attempt timed out after {:?}", t),
                    ))),
                },
                None => attempt_fut.await,
            };
            match result {
                Ok(resp) => return Ok(resp),
                Err(e) => {
                    let retryable = match &e {
                        PredictError::Grpc(s) => {
                            RetryPolicy::is_retryable_grpc(s.code())
                        }
                        PredictError::NonOk { status, .. } => {
                            RetryPolicy::is_retryable_status(*status)
                        }
                        _ => false,
                    };
                    if !retryable || attempt == self.retry_policy.max_retries {
                        return Err(e);
                    }
                    tracing::debug!(
                        "predict attempt {} failed (retryable): {}; retrying",
                        attempt,
                        e
                    );
                    last_err = Some(e);
                }
            }
        }

        Err(last_err.unwrap_or_else(|| {
            PredictError::Grpc(tonic::Status::internal("retry loop exited unexpectedly"))
        }))
    }

    async fn try_predict_once(
        &mut self,
        req: PredictRequest,
        expected_model_sha256: Option<&str>,
    ) -> Result<PredictResponse, PredictError> {
        let request_id = req.request_id.clone();
        let mut resp = self.inner.predict(req).await?.into_inner();
        if resp.status != PredictStatus::Ok as i32 {
            return Err(PredictError::NonOk {
                status: PredictStatus::try_from(resp.status)
                    .unwrap_or(PredictStatus::Unspecified),
                code: resp.error_code,
                message: resp.error_message,
            });
        }
        if let Some(expected) = expected_model_sha256 {
            if resp.model_sha256 != expected {
                return Err(PredictError::ModelShaMismatch {
                    expected: expected.to_string(),
                    got: resp.model_sha256.clone(),
                });
            }
        }

        let server_encrypted = !resp.y_aead_nonce.is_empty();
        match (self.aead_key.as_ref(), server_encrypted) {
            (Some(key), true) => {
                let aad = crate::aead::response_aad(
                    &request_id,
                    resp.y_dtype,
                    &resp.y_shape,
                    &resp.y_sha256,
                );
                let plaintext = crate::aead::open(
                    key,
                    &resp.y_payload,
                    &resp.y_aead_nonce,
                    &aad,
                )
                .map_err(|e| {
                    PredictError::Grpc(tonic::Status::data_loss(format!(
                        "aead open: {e}"
                    )))
                })?;
                resp.y_payload = plaintext.into();
                resp.y_aead_nonce = Default::default();
            }
            (Some(_), false) => {
                return Err(PredictError::Grpc(tonic::Status::failed_precondition(
                    "client expects AEAD but server response is plaintext",
                )));
            }
            (None, true) => {
                return Err(PredictError::Grpc(tonic::Status::failed_precondition(
                    "server returned AEAD ciphertext but client has no key",
                )));
            }
            (None, false) => {  }
        }
        let computed = sha256_hex(&resp.y_payload);
        if !resp.y_sha256.is_empty() && resp.y_sha256 != computed {
            return Err(PredictError::Sha256Mismatch {
                declared: resp.y_sha256.clone(),
                computed,
            });
        }
        Ok(resp)
    }

    fn next_backoff(&self, attempt: u32) -> Duration {
        self.retry_policy.next_backoff_for(attempt)
    }

    pub async fn health(&mut self) -> Result<crate::proto::HealthResponse, PredictError> {
        let resp = self.inner.health(HealthRequest {}).await?.into_inner();
        Ok(resp)
    }
}

async fn connect_endpoint_with_retry(
    ep: tonic::transport::Endpoint,
    retry_policy: RetryPolicy,
) -> Result<Channel, PredictError> {
    let mut last_err: Option<tonic::transport::Error> = None;
    for attempt in 0..=retry_policy.max_retries {
        if attempt > 0 {
            let delay = retry_policy.next_backoff_for(attempt - 1);
            tokio::time::sleep(delay).await;
        }
        match ep.connect().await {
            Ok(chan) => return Ok(chan),
            Err(e) => {
                tracing::debug!("connect attempt {} failed: {}", attempt, e);
                last_err = Some(e);
            }
        }
    }
    Err(PredictError::Transport(
        last_err.expect("at least one attempt"),
    ))
}

pub fn make_request_i32(
    request_id: &str,
    model_id: &str,
    shape: Vec<i64>,
    data: &[i32],
) -> PredictRequest {
    let payload = codec::encode_i32(data);
    let payload_sha256 = sha256_hex(&payload);
    PredictRequest {
        request_id: request_id.to_string(),
        model_id: model_id.to_string(),
        dtype: TensorDType::TensorDtypeInt32 as i32,
        shape,
        payload: payload.into(),
        payload_sha256,
        aead_nonce: Default::default(),
        payload_encoding: String::new(),
    }
}

pub fn make_request_f32(
    request_id: &str,
    model_id: &str,
    shape: Vec<i64>,
    data: &[f32],
) -> PredictRequest {
    let payload = codec::encode_f32(data);
    let payload_sha256 = sha256_hex(&payload);
    PredictRequest {
        request_id: request_id.to_string(),
        model_id: model_id.to_string(),
        dtype: TensorDType::TensorDtypeFloat32 as i32,
        shape,
        payload: payload.into(),
        payload_sha256,
        aead_nonce: Default::default(),
        payload_encoding: String::new(),
    }
}

pub fn make_request_i8(
    request_id: &str,
    model_id: &str,
    shape: Vec<i64>,
    data: &[i8],
) -> PredictRequest {
    let payload = codec::encode_i8(data);
    let payload_sha256 = sha256_hex(&payload);
    PredictRequest {
        request_id: request_id.to_string(),
        model_id: model_id.to_string(),
        dtype: TensorDType::TensorDtypeInt8 as i32,
        shape,
        payload: payload.into(),
        payload_sha256,
        aead_nonce: Default::default(),
        payload_encoding: String::new(),
    }
}
