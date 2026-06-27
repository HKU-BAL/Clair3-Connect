use anyhow::{Context, Result};
use clair3_connect_infer::device::Device;
use clair3_connect_infer::model_registry::Registry;
use clair3_connect_server::backpressure::BackpressureConfig;
use clair3_connect_server::service::PredictService;
use clair3_connect_transport::aead::{AeadKey, AEAD_KEY_BYTES};
use clair3_connect_transport::proto::clair3_connect_inference_server::Clair3ConnectInferenceServer;
use clair3_connect_transport::tls::ServerTlsOptions;
use std::path::PathBuf;
use std::time::Duration;
use tonic::transport::Server;

#[derive(Debug)]
struct Args {

    bind: Option<String>,

    registry: Option<PathBuf>,

    devices: Option<Vec<Device>>,
    max_inflight: usize,

    cpu_workers: usize,

    max_batch_size: usize,

    flush_timeout_ms: u64,
    retry_after_ms: u64,
    tls_cert: Option<PathBuf>,
    tls_key: Option<PathBuf>,
    tls_client_ca: Option<PathBuf>,

    metrics_bind: Option<String>,

    aead_key_path: Option<PathBuf>,

    force_interactive: bool,

    force_non_interactive: bool,

    pid_file: Option<PathBuf>,
}

fn parse_args() -> Args {
    let mut bind: Option<String> = None;
    let mut registry: Option<PathBuf> = None;
    let mut device_spec: Option<String> = None;
    let mut max_inflight: usize = BackpressureConfig::default().max_inflight;
    let mut cpu_workers: usize = 1;
    let bc_default = clair3_connect_infer::batcher::BatchConfig::default();
    let mut max_batch_size: usize = bc_default.max_batch_size;
    let mut flush_timeout_ms: u64 = bc_default.flush_timeout.as_millis() as u64;
    let mut retry_after_ms: u64 = BackpressureConfig::default().retry_after.as_millis() as u64;
    let mut tls_cert: Option<PathBuf> = None;
    let mut tls_key: Option<PathBuf> = None;
    let mut tls_client_ca: Option<PathBuf> = None;
    let mut force_interactive = false;
    let mut force_non_interactive = false;
    let mut metrics_bind: Option<String> = None;
    let mut aead_key_path: Option<PathBuf> = None;
    let mut pid_file: Option<PathBuf> = None;

    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--bind" => bind = Some(args.next().expect("--bind needs a value")),
            "--model-registry" => {
                registry = Some(PathBuf::from(
                    args.next().expect("--model-registry needs a value"),
                ))
            }
            "--device" => device_spec = Some(args.next().expect("--device needs a value")),
            "--interactive" | "-i" => force_interactive = true,
            "--non-interactive" => force_non_interactive = true,
            "--max-inflight" => {
                max_inflight = args
                    .next()
                    .expect("--max-inflight needs a value")
                    .parse()
                    .expect("--max-inflight must be a non-negative integer")
            }
            "--cpu-workers" => {
                cpu_workers = args
                    .next()
                    .expect("--cpu-workers needs a value")
                    .parse()
                    .expect("--cpu-workers must be a positive integer");
                assert!(cpu_workers >= 1, "--cpu-workers must be >= 1");
            }
            "--max-batch-size" => {
                max_batch_size = args
                    .next()
                    .expect("--max-batch-size needs a value")
                    .parse()
                    .expect("--max-batch-size must be a positive integer");
                assert!(max_batch_size >= 1, "--max-batch-size must be >= 1");
            }
            "--flush-timeout-ms" => {
                flush_timeout_ms = args
                    .next()
                    .expect("--flush-timeout-ms needs a value")
                    .parse()
                    .expect("--flush-timeout-ms must be u64 ms")
            }
            "--retry-after-ms" => {
                retry_after_ms = args
                    .next()
                    .expect("--retry-after-ms needs a value")
                    .parse()
                    .expect("--retry-after-ms must be u64 ms")
            }
            "--tls-cert" => {
                tls_cert = Some(PathBuf::from(args.next().expect("--tls-cert PATH")))
            }
            "--tls-key" => {
                tls_key = Some(PathBuf::from(args.next().expect("--tls-key PATH")))
            }
            "--tls-client-ca" => {
                tls_client_ca = Some(PathBuf::from(
                    args.next().expect("--tls-client-ca PATH"),
                ))
            }
            "--metrics-bind" => {
                metrics_bind = Some(args.next().expect("--metrics-bind ADDR"))
            }
            "--aead-key" => {
                aead_key_path = Some(PathBuf::from(args.next().expect("--aead-key PATH")))
            }
            "--pid-file" => {
                pid_file = Some(PathBuf::from(args.next().expect("--pid-file PATH")))
            }
            "-h" | "--help" => {
                eprintln!(
                    "clair3-connect-server [--bind ADDR] [--model-registry PATH] \
                     [--device auto|cpu|cuda:N|cuda:N,cuda:M,...] [--max-inflight N] [--cpu-workers N] \
                     [--max-batch-size N (64)] [--flush-timeout-ms MS (4)] [--retry-after-ms MS] \
                     [--tls-cert PEM --tls-key PEM [--tls-client-ca PEM]] \
                     [--metrics-bind ADDR] [--aead-key PATH (32-byte binary)] \
                     [--pid-file PATH] [--interactive | --non-interactive]

Live control panel (reconfigure a RUNNING server, no restart):
  clair3-connect-server control --model-registry PATH [--pid-file PATH | --pid N]
      Opens the add/edit/remove model menu, then SIGHUPs the server to
      hot-reload its models. Start the server with --pid-file so this can
      find it. (You can also just `kill -HUP <pid>` after editing the YAML.)

Defaults (when not specified):
  --bind            127.0.0.1:50051
  --model-registry  models/model_registry.example.yaml
  --device          AUTO — all CUDA GPUs if any detected, else cpu
  --max-inflight    8 (per device)
  --retry-after-ms  100

Tip: run with no flags (on a TTY) to start a guided setup that walks
you through bind, registry, GPUs, backpressure, and metrics, then
security. TLS + mTLS and payload AEAD default to ON — the wizard can
mint a self-signed PKI (CA + server + client) and a 32-byte AEAD key
on the spot, so a secured server is one Enter away (opt out
explicitly). The registry step opens a model-management menu: add /
edit / remove / verify entries, saved back to YAML.

Note: --tls-cert/--tls-key/--aead-key and the non-interactive path are
unchanged and stay opt-in; the default-on behaviour is the wizard's.

The binary embeds a libtorch search path (DT_RUNPATH); no need to
export LD_LIBRARY_PATH unless you want to override the default."
                );
                std::process::exit(0);
            }
            other => {
                eprintln!("unknown flag: {other}");
                std::process::exit(1);
            }
        }
    }
    if tls_cert.is_some() != tls_key.is_some() {
        eprintln!("--tls-cert and --tls-key must be specified together");
        std::process::exit(1);
    }

    if tls_client_ca.is_some() && (tls_cert.is_none() || tls_key.is_none()) {
        eprintln!(
            "--tls-client-ca requires --tls-cert + --tls-key (mTLS cannot run without a server identity)"
        );
        std::process::exit(1);
    }
    let devices: Option<Vec<Device>> = device_spec.as_deref().map(|s| {
        Device::parse_list(s).unwrap_or_else(|e| panic!("bad --device {s:?}: {e}"))
    });
    Args {
        bind,
        registry,
        devices,
        max_inflight,
        cpu_workers,
        max_batch_size,
        flush_timeout_ms,
        retry_after_ms,
        tls_cert,
        tls_key,
        tls_client_ca,
        metrics_bind,
        aead_key_path,
        force_interactive,
        force_non_interactive,
        pid_file,
    }
}

fn auto_detect_devices() -> Vec<Device> {
    let _ = clair3_connect_infer::_ensure_cuda_linked();
    if tch::Cuda::is_available() {
        let n = tch::Cuda::device_count();
        if n > 0 {
            tracing::info!(
                n_gpus = n,
                "auto-selected all CUDA devices (set --device to override)"
            );
            return (0..n as usize).map(Device::Cuda).collect();
        }
    }
    tracing::info!("no CUDA detected — auto-selected CPU (set --device to override)");
    vec![Device::Cpu]
}

fn should_run_interactive(args: &Args) -> bool {
    if args.force_non_interactive {
        return false;
    }
    if args.force_interactive {
        return true;
    }

    let nothing_passed = std::env::args().len() <= 1;
    nothing_passed && console::user_attended()
}

fn merge_interactive(
    args: &mut Args,
    iv: clair3_connect_server::interactive::InteractiveServerArgs,
) {
    if args.bind.is_none() {
        args.bind = Some(iv.bind);
    }
    if args.registry.is_none() {
        args.registry = Some(iv.registry);
    }
    if args.devices.is_none() {
        args.devices = Some(iv.devices);
    }

    if args.max_inflight == BackpressureConfig::default().max_inflight {
        args.max_inflight = iv.max_inflight;
    }
    if args.retry_after_ms == BackpressureConfig::default().retry_after.as_millis() as u64 {
        args.retry_after_ms = iv.retry_after_ms;
    }
    if args.metrics_bind.is_none() {
        args.metrics_bind = iv.metrics_bind;
    }
    if args.tls_cert.is_none() {
        args.tls_cert = iv.tls_cert;
    }
    if args.tls_key.is_none() {
        args.tls_key = iv.tls_key;
    }
    if args.tls_client_ca.is_none() {
        args.tls_client_ca = iv.tls_client_ca;
    }
    if args.aead_key_path.is_none() {
        args.aead_key_path = iv.aead_key_path;
    }
}

fn run_control(args: &[String]) -> Result<()> {
    let mut registry: Option<PathBuf> = None;
    let mut pid_file: Option<PathBuf> = None;
    let mut pid: Option<i32> = None;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--model-registry" => registry = it.next().map(PathBuf::from),
            "--pid-file" => pid_file = it.next().map(PathBuf::from),
            "--pid" => pid = it.next().and_then(|s| s.parse().ok()),
            "-h" | "--help" => {
                eprintln!(
                    "clair3-connect-server control --model-registry PATH [--pid-file PATH | --pid N]\n\n\
                     Opens the add/edit/remove model menu against the registry, then sends\n\
                     SIGHUP to the running server so it hot-reloads its models (no restart).\n\
                     Point --model-registry at the SAME yaml the server was started with;\n\
                     start that server with --pid-file so this command can signal it."
                );
                return Ok(());
            }
            other => anyhow::bail!("unknown `control` argument: {other} (try: control --help)"),
        }
    }
    let registry = registry.context("control: --model-registry <yaml> is required")?;
    let theme = dialoguer::theme::ColorfulTheme::default();
    clair3_connect_server::registry_admin::manage(&theme, &registry)?;

    let pid = pid.or_else(|| {
        pid_file
            .as_ref()
            .and_then(|p| std::fs::read_to_string(p).ok())
            .and_then(|s| s.trim().parse().ok())
    });
    match pid {
        Some(p) => {
            let ok = std::process::Command::new("kill")
                .args(["-HUP", &p.to_string()])
                .status()
                .map(|s| s.success())
                .unwrap_or(false);
            if ok {
                println!("✓ Sent SIGHUP to server pid {p} — model registry reloaded live.");
            } else {
                println!("Registry edited, but signaling pid {p} failed. Run `kill -HUP {p}` manually.");
            }
        }
        None => println!(
            "Registry edited. To apply on the running server, send it SIGHUP:\n  \
             kill -HUP <server-pid>\n(or re-run with --pid-file/--pid to signal it automatically)."
        ),
    }
    Ok(())
}

#[tokio::main]
async fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();

    if std::env::var_os("PYTORCH_TENSOREXPR").is_none() {
        std::env::set_var("PYTORCH_TENSOREXPR", "0");
        tracing::debug!("PYTORCH_TENSOREXPR=0 set (disables the NVRTC fuser; export =1 to re-enable)");
    }

    let raw: Vec<String> = std::env::args().collect();
    if raw.get(1).map(String::as_str) == Some("control") {
        return run_control(&raw[2..]);
    }

    let mut args = parse_args();

    if should_run_interactive(&args) {
        let defaults = clair3_connect_server::interactive::InteractiveServerDefaults {
            bind: args.bind.clone(),
            registry: args.registry.clone(),
            device_spec: args
                .devices
                .as_ref()
                .map(|v| v.iter().map(|d| d.to_string()).collect::<Vec<_>>().join(",")),
            max_inflight: Some(args.max_inflight),
            retry_after_ms: Some(args.retry_after_ms),
            metrics_bind: args.metrics_bind.clone(),
            tls_cert: args.tls_cert.clone(),
            tls_key: args.tls_key.clone(),
            tls_client_ca: args.tls_client_ca.clone(),
            aead_key_path: args.aead_key_path.clone(),
        };
        let iv = clair3_connect_server::interactive::run(defaults)
            .context("interactive setup")?;
        merge_interactive(&mut args, iv);
    }

    let bind = args
        .bind
        .clone()
        .unwrap_or_else(|| "127.0.0.1:50051".to_string());
    let registry_path = args.registry.clone().unwrap_or_else(|| {
        PathBuf::from("models/model_registry.example.yaml")
    });
    let devices = match args.devices.clone() {
        Some(d) => d,
        None => auto_detect_devices(),
    };

    let devices = if args.cpu_workers > 1 && devices.iter().all(|d| matches!(d, Device::Cpu)) {
        let cores = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(8);

        let per = std::env::var("OMP_NUM_THREADS")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or_else(|| (cores * 3 / 4 / args.cpu_workers).max(1));
        std::env::set_var("OMP_NUM_THREADS", per.to_string());
        tch::set_num_threads(per as i32);
        tracing::info!(
            cpu_workers = args.cpu_workers,
            threads_per_worker = per,
            "parallel CPU inference: {} cpu batchers",
            args.cpu_workers
        );
        std::iter::repeat(Device::Cpu)
            .take(args.cpu_workers)
            .collect()
    } else {
        devices
    };

    tracing::info!(bind = %bind, registry = %registry_path.display(), "loading model registry");
    let reg = Registry::load(&registry_path).context("load registry")?;
    let yaml_dir = registry_path.parent().unwrap().to_path_buf();
    let bp = BackpressureConfig {
        max_inflight: args.max_inflight,
        retry_after: Duration::from_millis(args.retry_after_ms),
    };

    let tch_devices: Vec<tch::Device> = devices.iter().map(|d| d.as_tch()).collect();
    let device_count = tch_devices.len();
    let mut svc = PredictService::load_with_options(
        &reg,
        &yaml_dir,
        &tch_devices,
        bp,
        clair3_connect_infer::batcher::BatchConfig {
            max_batch_size: args.max_batch_size,
            flush_timeout: Duration::from_millis(args.flush_timeout_ms),
            ..clair3_connect_infer::batcher::BatchConfig::default()
        },
    )
    .context("load models")?;
    tracing::info!(
        device_count,
        "PredictService loaded ({} devices: {:?})",
        device_count,
        devices
    );
    tracing::info!(
        max_inflight = svc.max_inflight(),
        retry_after_ms = args.retry_after_ms,
        "backpressure gate configured"
    );
    tracing::info!(
        max_batch_size = args.max_batch_size,
        flush_timeout_ms = args.flush_timeout_ms,
        "dynamic batcher configured"
    );

    if let Some(p) = args.aead_key_path.as_ref() {
        let bytes = std::fs::read(p).with_context(|| format!("read {}", p.display()))?;
        if bytes.len() != AEAD_KEY_BYTES {
            anyhow::bail!(
                "--aead-key file {} must be exactly {} bytes (got {})",
                p.display(),
                AEAD_KEY_BYTES,
                bytes.len()
            );
        }
        let key = AeadKey::from_bytes(&bytes).context("AeadKey::from_bytes")?;
        svc.set_aead_key(key);
        tracing::info!(path = %p.display(), "AEAD enabled (AES-256-GCM)");
    }

    let addr = bind.parse().context("parse --bind addr")?;

    if let Some(mb) = args.metrics_bind.as_deref() {
        let metrics_addr = mb.parse().context("parse --metrics-bind addr")?;
        clair3_connect_server::metrics_server::spawn(metrics_addr, svc.clone())
            .await
            .context("metrics_server::spawn")?;
    }

    if let Some(pf) = args.pid_file.as_ref() {
        std::fs::write(pf, std::process::id().to_string())
            .with_context(|| format!("write pid file {}", pf.display()))?;
        tracing::info!(pid_file = %pf.display(), pid = std::process::id(), "wrote pid file");
    }

    {
        let svc = svc.clone();
        let reg = registry_path.clone();
        tokio::spawn(async move {
            let mut hup =
                match tokio::signal::unix::signal(tokio::signal::unix::SignalKind::hangup()) {
                    Ok(s) => s,
                    Err(e) => {
                        tracing::error!("SIGHUP handler unavailable, live reload disabled: {e}");
                        return;
                    }
                };
            while hup.recv().await.is_some() {
                tracing::info!(registry = %reg.display(), "SIGHUP — reloading model registry");
                match svc.reload_from_registry_file(&reg) {
                    Ok(s) => tracing::info!(
                        added = ?s.added, updated = ?s.updated,
                        removed = ?s.removed, unchanged = s.unchanged,
                        "model registry reloaded live"
                    ),
                    Err(e) => {
                        tracing::error!("registry reload failed (keeping current models): {e:#}")
                    }
                }
            }
        });
    }

    let mut builder = Server::builder().tcp_nodelay(true);
    if let (Some(cert), Some(key)) = (args.tls_cert.as_ref(), args.tls_key.as_ref()) {
        let opts = ServerTlsOptions::from_paths(cert, key, args.tls_client_ca.as_deref())
            .context("load server TLS")?;
        let mutual = opts.is_mutual();
        builder = builder
            .tls_config(opts.build())
            .context("tonic tls_config")?;
        tracing::info!(
            mutual,
            cert = %cert.display(),
            "TLS enabled"
        );
    } else {
        tracing::warn!("starting in plaintext mode — TLS NOT enabled");
    }
    tracing::info!("clair3-connect-server listening on {}", addr);

    builder
        .add_service(
            Clair3ConnectInferenceServer::new(svc)
                .accept_compressed(tonic::codec::CompressionEncoding::Zstd)
                .accept_compressed(tonic::codec::CompressionEncoding::Gzip)
                .send_compressed(tonic::codec::CompressionEncoding::Zstd),
        )
        .serve(addr)
        .await?;
    Ok(())
}
