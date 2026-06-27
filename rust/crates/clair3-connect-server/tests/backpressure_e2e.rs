use clair3_connect_infer::model_registry::Registry;
use clair3_connect_server::backpressure::BackpressureConfig;
use clair3_connect_server::service::PredictService;
use clair3_connect_transport::predict_client::{make_request_i32, PredictClient, PredictError, RetryPolicy};
use clair3_connect_transport::proto::clair3_connect_inference_server::Clair3ConnectInferenceServer;
use std::path::PathBuf;
use std::time::Duration;
use tokio::net::TcpListener;
use tonic::transport::Server;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

#[tokio::test]
async fn server_rejects_when_gate_is_zero() {
    let yaml = workspace_root().join("models/model_registry.example.yaml");
    if !yaml.exists() {
        eprintln!("Skipping — model registry yaml missing");
        return;
    }
    let reg = Registry::load(&yaml).expect("registry");
    let yaml_dir = yaml.parent().unwrap().to_path_buf();

    let bp = BackpressureConfig {
        max_inflight: 0,
        retry_after: Duration::from_millis(250),
    };
    let svc = match PredictService::load_with_backpressure(&reg, &yaml_dir, tch::Device::Cpu, bp) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Skipping — service load failed (TorchScript fixture?): {e}");
            return;
        }
    };
    assert_eq!(svc.max_inflight(), 0);

    let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    let incoming = tokio_stream::wrappers::TcpListenerStream::new(listener);
    let server = tokio::spawn(async move {
        Server::builder()
            .add_service(Clair3ConnectInferenceServer::new(svc))
            .serve_with_incoming(incoming)
            .await
            .ok();
    });

    let endpoint = format!("http://{}", addr);
    let mut client = PredictClient::from_channel(
        tonic::transport::Endpoint::from_shared(endpoint)
            .unwrap()
            .connect()
            .await
            .expect("connect"),
        RetryPolicy::none(),
    );

    let req = make_request_i32("bp1", "irrelevant", vec![1, 33, 18], &vec![0_i32; 33 * 18]);
    let err = client.predict(req).await.expect_err("must be rejected");
    match err {
        PredictError::Grpc(status) => {
            assert_eq!(status.code(), tonic::Code::ResourceExhausted, "code");
            let hint = status
                .metadata()
                .get("retry-after-ms")
                .expect("retry-after-ms header present")
                .to_str()
                .unwrap()
                .to_string();
            assert_eq!(hint, "250", "retry-after-ms value");
        }
        other => panic!("expected Grpc(ResourceExhausted), got {:?}", other),
    }

    server.abort();
}
