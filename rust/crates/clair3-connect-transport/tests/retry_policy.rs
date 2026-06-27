use clair3_connect_transport::predict_client::{
    make_request_i8, PredictClient, PredictError, RetryPolicy,
};
use clair3_connect_transport::proto::clair3_connect_inference_server::{
    Clair3ConnectInference, Clair3ConnectInferenceServer,
};
use clair3_connect_transport::proto::{
    HealthRequest, HealthResponse, PredictRequest, PredictResponse, PredictStatus,
};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::Duration;
use tokio::net::TcpListener;
use tonic::transport::{Channel, Server};
use tonic::{Request, Response, Status};

#[derive(Debug, Clone)]
enum Reply {

    GrpcError(tonic::Code),

    Typed(PredictStatus),

    Ok { y_payload: Vec<u8> },

    Stall(Duration),
}

#[derive(Clone)]
struct MockPredictService {

    sequence: Arc<Vec<Reply>>,

    call_count: Arc<AtomicUsize>,
}

impl MockPredictService {
    fn new(sequence: Vec<Reply>) -> Self {
        Self {
            sequence: Arc::new(sequence),
            call_count: Arc::new(AtomicUsize::new(0)),
        }
    }

    #[allow(dead_code)]
    fn calls(&self) -> usize {
        self.call_count.load(Ordering::SeqCst)
    }
}

#[tonic::async_trait]
impl Clair3ConnectInference for MockPredictService {
    async fn predict(
        &self,
        request: Request<PredictRequest>,
    ) -> Result<Response<PredictResponse>, Status> {
        let n = self.call_count.fetch_add(1, Ordering::SeqCst);
        let idx = n.min(self.sequence.len() - 1);
        let reply = &self.sequence[idx];
        let req = request.into_inner();
        match reply {
            Reply::GrpcError(code) => Err(Status::new(*code, format!("mock error: {n}"))),
            Reply::Typed(status) => Ok(Response::new(PredictResponse {
                request_id: req.request_id,
                status: *status as i32,
                error_code: "mock_retry".to_string(),
                error_message: format!("mock typed status at call {n}"),
                ..Default::default()
            })),
            Reply::Ok { y_payload } => Ok(Response::new(PredictResponse {
                request_id: req.request_id,
                status: PredictStatus::Ok as i32,
                error_code: String::new(),
                error_message: String::new(),
                y_dtype: 0,
                y_shape: vec![],
                y_payload: y_payload.clone().into(),
                y_sha256: clair3_connect_transport::checksum::sha256_hex(y_payload),
                input_payload_sha256: req.payload_sha256,
                model_sha256: "mock_model_sha".to_string(),
                y_aead_nonce: Default::default(),
            })),
            Reply::Stall(d) => {
                tokio::time::sleep(*d).await;
                Ok(Response::new(PredictResponse {
                    request_id: req.request_id,
                    status: PredictStatus::Ok as i32,
                    ..Default::default()
                }))
            }
        }
    }

    async fn health(
        &self,
        _request: Request<HealthRequest>,
    ) -> Result<Response<HealthResponse>, Status> {
        Ok(Response::new(HealthResponse {
            server_version: "mock".to_string(),
            registered_models: Default::default(),
        }))
    }
}

async fn spawn_mock_server(svc: MockPredictService) -> (Channel, tokio::task::JoinHandle<()>) {
    let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    let incoming = tokio_stream::wrappers::TcpListenerStream::new(listener);
    let handle = tokio::spawn(async move {
        Server::builder()
            .add_service(Clair3ConnectInferenceServer::new(svc))
            .serve_with_incoming(incoming)
            .await
            .ok();
    });
    tokio::task::yield_now().await;
    let channel = Channel::from_shared(format!("http://{}", addr))
        .unwrap()
        .connect()
        .await
        .unwrap();
    (channel, handle)
}

#[tokio::test]
async fn succeeds_after_two_unavailable_retries() {
    let svc = MockPredictService::new(vec![
        Reply::GrpcError(tonic::Code::Unavailable),
        Reply::GrpcError(tonic::Code::Unavailable),
        Reply::Ok { y_payload: vec![1, 2, 3, 4] },
    ]);
    let calls_ref = svc.call_count.clone();
    let (chan, server) = spawn_mock_server(svc).await;

    let mut client = PredictClient::from_channel(
        chan,
        RetryPolicy {
            max_retries: 3,
            initial_delay: Duration::from_millis(10),
            max_delay: Duration::from_millis(100),
            per_attempt_timeout: None,
        },
    );
    let req = make_request_i8("r1", "m1", vec![1, 4], &[1, 2, 3, 4]);
    let resp = client.predict(req).await.expect("should succeed after retries");
    assert_eq!(resp.y_payload, vec![1, 2, 3, 4]);
    assert_eq!(calls_ref.load(Ordering::SeqCst), 3, "expected 3 attempts");

    server.abort();
}

#[tokio::test]
async fn typed_predict_status_retry_is_retried() {
    let svc = MockPredictService::new(vec![
        Reply::Typed(PredictStatus::Retry),
        Reply::Ok { y_payload: vec![0xAB] },
    ]);
    let calls_ref = svc.call_count.clone();
    let (chan, server) = spawn_mock_server(svc).await;
    let mut client = PredictClient::from_channel(
        chan,
        RetryPolicy {
            max_retries: 2,
            initial_delay: Duration::from_millis(10),
            max_delay: Duration::from_millis(50),
            per_attempt_timeout: None,
        },
    );
    let req = make_request_i8("r2", "m1", vec![1, 1], &[1]);
    let resp = client.predict(req).await.expect("should succeed after typed Retry");
    assert_eq!(resp.y_payload, vec![0xAB]);
    assert_eq!(calls_ref.load(Ordering::SeqCst), 2);
    server.abort();
}

#[tokio::test]
async fn internal_error_is_fatal_no_retry() {
    let svc = MockPredictService::new(vec![
        Reply::GrpcError(tonic::Code::Internal),
        Reply::Ok { y_payload: vec![1] },
    ]);
    let calls_ref = svc.call_count.clone();
    let (chan, server) = spawn_mock_server(svc).await;
    let mut client = PredictClient::from_channel(
        chan,
        RetryPolicy {
            max_retries: 3,
            initial_delay: Duration::from_millis(10),
            max_delay: Duration::from_millis(50),
            per_attempt_timeout: None,
        },
    );
    let req = make_request_i8("r3", "m1", vec![1, 1], &[1]);
    let result = client.predict(req).await;
    assert!(matches!(result, Err(PredictError::Grpc(_))));
    assert_eq!(calls_ref.load(Ordering::SeqCst), 1, "Internal must not retry");
    server.abort();
}

#[tokio::test]
async fn exhausts_max_retries_then_returns_last_error() {

    let svc = MockPredictService::new(vec![Reply::GrpcError(tonic::Code::Unavailable)]);
    let calls_ref = svc.call_count.clone();
    let (chan, server) = spawn_mock_server(svc).await;
    let mut client = PredictClient::from_channel(
        chan,
        RetryPolicy {
            max_retries: 2,
            initial_delay: Duration::from_millis(5),
            max_delay: Duration::from_millis(20),
            per_attempt_timeout: None,
        },
    );
    let req = make_request_i8("r4", "m1", vec![1, 1], &[1]);
    let result = client.predict(req).await;
    assert!(matches!(result, Err(PredictError::Grpc(_))));

    assert_eq!(calls_ref.load(Ordering::SeqCst), 3);
    server.abort();
}

#[tokio::test]
async fn aborted_and_resource_exhausted_retry() {

    let svc = MockPredictService::new(vec![
        Reply::GrpcError(tonic::Code::Aborted),
        Reply::GrpcError(tonic::Code::ResourceExhausted),
        Reply::Ok { y_payload: vec![0xCC] },
    ]);
    let calls_ref = svc.call_count.clone();
    let (chan, server) = spawn_mock_server(svc).await;
    let mut client = PredictClient::from_channel(
        chan,
        RetryPolicy {
            max_retries: 5,
            initial_delay: Duration::from_millis(5),
            max_delay: Duration::from_millis(30),
            per_attempt_timeout: None,
        },
    );
    let req = make_request_i8("r5", "m1", vec![1, 1], &[1]);
    let resp = client.predict(req).await.expect("should succeed");
    assert_eq!(resp.y_payload, vec![0xCC]);
    assert_eq!(calls_ref.load(Ordering::SeqCst), 3);
    server.abort();
}

#[test]
fn retry_policy_base_delay_doubles() {
    let p = RetryPolicy {
        max_retries: 10,
        initial_delay: Duration::from_millis(100),
        max_delay: Duration::from_secs(60),
        per_attempt_timeout: None,
    };
    assert_eq!(p.base_delay_for(0), Duration::from_millis(100));
    assert_eq!(p.base_delay_for(1), Duration::from_millis(200));
    assert_eq!(p.base_delay_for(2), Duration::from_millis(400));
    assert_eq!(p.base_delay_for(3), Duration::from_millis(800));
}

#[test]
fn retry_policy_caps_at_max_delay() {
    let p = RetryPolicy {
        max_retries: 20,
        initial_delay: Duration::from_millis(100),
        max_delay: Duration::from_secs(1),
        per_attempt_timeout: None,
    };

    assert_eq!(p.base_delay_for(4), Duration::from_secs(1));
    assert_eq!(p.base_delay_for(15), Duration::from_secs(1));
}

#[tokio::test]
async fn per_attempt_timeout_triggers_retry() {

    let svc = MockPredictService::new(vec![Reply::Stall(Duration::from_millis(200))]);
    let calls_ref = svc.call_count.clone();
    let (chan, server) = spawn_mock_server(svc).await;
    let mut client = PredictClient::from_channel(
        chan,
        RetryPolicy {
            max_retries: 2,
            initial_delay: Duration::from_millis(5),
            max_delay: Duration::from_millis(20),
            per_attempt_timeout: Some(Duration::from_millis(50)),
        },
    );
    let req = make_request_i8("r6", "m1", vec![1, 1], &[1]);
    let result = client.predict(req).await;

    assert!(matches!(result, Err(PredictError::Grpc(_))));
    assert_eq!(calls_ref.load(Ordering::SeqCst), 1);
    server.abort();
}

#[tokio::test]
async fn connect_with_retry_retries_on_refused() {

    use std::time::Instant;
    let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    drop(listener);
    let endpoint = format!("http://{}", addr);

    let policy = RetryPolicy {
        max_retries: 2,

        initial_delay: Duration::from_millis(30),
        max_delay: Duration::from_millis(200),
        per_attempt_timeout: None,
    };
    let start = Instant::now();
    let res = PredictClient::connect_with_retry(endpoint, policy).await;
    let elapsed = start.elapsed();

    assert!(
        elapsed >= Duration::from_millis(50),
        "elapsed {:?} too short — retry loop not running?",
        elapsed
    );

    let err = res.err().expect("connect must fail after retries exhausted");
    assert!(
        matches!(err, PredictError::Transport(_)),
        "expected PredictError::Transport, got {}",
        err
    );
}

#[tokio::test]
async fn connect_with_retry_fast_path_when_max_retries_zero() {

    use std::time::Instant;
    let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    drop(listener);
    let endpoint = format!("http://{}", addr);

    let policy = RetryPolicy {
        max_retries: 0,
        initial_delay: Duration::from_secs(10),
        max_delay: Duration::from_secs(20),
        per_attempt_timeout: None,
    };
    let start = Instant::now();
    let _ = PredictClient::connect_with_retry(endpoint, policy).await;
    let elapsed = start.elapsed();
    assert!(
        elapsed < Duration::from_secs(1),
        "fast-path should not sleep: elapsed {:?}",
        elapsed
    );
}

#[test]
fn retryable_classification() {
    assert!(RetryPolicy::is_retryable_grpc(tonic::Code::Unavailable));
    assert!(RetryPolicy::is_retryable_grpc(tonic::Code::Aborted));
    assert!(RetryPolicy::is_retryable_grpc(tonic::Code::ResourceExhausted));
    assert!(!RetryPolicy::is_retryable_grpc(tonic::Code::Internal));
    assert!(!RetryPolicy::is_retryable_grpc(tonic::Code::InvalidArgument));
    assert!(!RetryPolicy::is_retryable_grpc(tonic::Code::NotFound));

    assert!(RetryPolicy::is_retryable_status(PredictStatus::Retry));
    assert!(!RetryPolicy::is_retryable_status(PredictStatus::Failed));
    assert!(!RetryPolicy::is_retryable_status(PredictStatus::Rejected));
}
