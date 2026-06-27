use anyhow::{bail, Context, Result};
use console::{style, Term};
use dialoguer::{theme::ColorfulTheme, Confirm, Input, Select};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct InteractiveServerArgs {

    pub bind: String,

    pub registry: PathBuf,

    pub devices: Vec<clair3_connect_infer::device::Device>,

    pub max_inflight: usize,

    pub retry_after_ms: u64,

    pub metrics_bind: Option<String>,

    pub tls_cert: Option<PathBuf>,

    pub tls_key: Option<PathBuf>,

    pub tls_client_ca: Option<PathBuf>,

    pub aead_key_path: Option<PathBuf>,
}

#[allow(missing_docs)]
#[derive(Debug, Clone, Default)]
pub struct InteractiveServerDefaults {
    pub bind: Option<String>,
    pub registry: Option<PathBuf>,
    pub device_spec: Option<String>,
    pub max_inflight: Option<usize>,
    pub retry_after_ms: Option<u64>,
    pub metrics_bind: Option<String>,
    pub tls_cert: Option<PathBuf>,
    pub tls_key: Option<PathBuf>,
    pub tls_client_ca: Option<PathBuf>,
    pub aead_key_path: Option<PathBuf>,
}

pub fn run(defaults: InteractiveServerDefaults) -> Result<InteractiveServerArgs> {
    let theme = ColorfulTheme::default();
    let term = Term::stdout();

    print_banner(&term)?;

    let bind = step_text(
        &theme,
        1,
        "gRPC bind address",
        defaults.bind.as_deref().unwrap_or("127.0.0.1:50051"),
    )?;

    let registry = step_registry(
        &theme,
        2,
        defaults.registry.as_deref().unwrap_or(Path::new(
            "models/model_registry.example.yaml",
        )),
    )?;

    let devices = step_devices(&theme, 3, defaults.device_spec.as_deref())?;

    let (max_inflight, retry_after_ms) = step_backpressure(
        &theme,
        4,
        defaults.max_inflight.unwrap_or(8),
        defaults.retry_after_ms.unwrap_or(100),
        devices.len(),
    )?;

    let metrics_bind = step_metrics(&theme, 5, defaults.metrics_bind.as_deref())?;

    let (tls_cert, tls_key, tls_client_ca) = step_tls(
        &theme,
        6,
        &bind,
        defaults.tls_cert.as_deref(),
        defaults.tls_key.as_deref(),
        defaults.tls_client_ca.as_deref(),
    )?;

    let aead_key_path = step_aead(&theme, 7, defaults.aead_key_path.as_deref())?;

    let args = InteractiveServerArgs {
        bind,
        registry,
        devices,
        max_inflight,
        retry_after_ms,
        metrics_bind,
        tls_cert,
        tls_key,
        tls_client_ca,
        aead_key_path,
    };

    print_summary(&term, &args)?;
    if !Confirm::with_theme(&theme)
        .with_prompt("Start the server with these settings?")
        .default(true)
        .interact()?
    {
        bail!("setup cancelled by user");
    }
    Ok(args)
}

fn print_banner(term: &Term) -> Result<()> {
    let line = "─".repeat(56);
    term.write_line("")?;
    term.write_line(&format!("{}", style(&line).cyan()))?;
    term.write_line(&format!(
        "{}",
        style("  Clair3-connect Server — interactive setup").bold().cyan()
    ))?;
    term.write_line(&format!(
        "{}",
        style("  Re-run with --non-interactive to skip prompts").dim()
    ))?;
    term.write_line(&format!("{}", style(&line).cyan()))?;
    term.write_line("")?;
    Ok(())
}

fn step_header(step: u8, total: u8, prompt: &str) -> String {
    format!(
        "{} {}",
        style(format!("[Step {step}/{total}]")).bold().green(),
        style(prompt).bold()
    )
}

fn step_text(
    theme: &ColorfulTheme,
    step: u8,
    prompt: &str,
    default: &str,
) -> Result<String> {
    println!();
    println!("{}", step_header(step, 7, prompt));
    let value: String = Input::with_theme(theme)
        .with_prompt("  value")
        .default(default.to_string())
        .interact_text()
        .with_context(|| format!("read {prompt}"))?;
    Ok(value)
}

fn step_registry(
    theme: &ColorfulTheme,
    step: u8,
    default: &Path,
) -> Result<PathBuf> {
    println!();
    println!("{}", step_header(step, 7, "Model registry YAML"));
    let path = loop {
        let input: String = Input::with_theme(theme)
            .with_prompt("  path")
            .default(default.display().to_string())
            .interact_text()
            .context("read registry path")?;
        let p = PathBuf::from(shellexpand_tilde(&input));
        if !p.exists() {
            println!(
                "  {}: file does not exist — try again or Ctrl-C to abort",
                style("error").red()
            );
            continue;
        }
        break p;
    };

    crate::registry_admin::print_registry_status(&path);

    let manage = Confirm::with_theme(theme)
        .with_prompt("  Open the model management menu (add / edit / remove)?")
        .default(false)
        .interact()
        .context("confirm registry mgmt")?;
    if manage {
        crate::registry_admin::manage(theme, &path)?;
    }
    Ok(path)
}

fn step_devices(
    theme: &ColorfulTheme,
    step: u8,
    default_spec: Option<&str>,
) -> Result<Vec<clair3_connect_infer::device::Device>> {
    use clair3_connect_infer::device::Device;
    println!();
    println!("{}", step_header(step, 7, "Compute devices"));

    let _ = clair3_connect_infer::_ensure_cuda_linked();
    let cuda_count = if tch::Cuda::is_available() {
        tch::Cuda::device_count()
    } else {
        0
    };
    if cuda_count > 0 {
        println!(
            "  {} {}",
            style("CUDA detected:").green(),
            style(format!("{cuda_count} device(s)")).bold()
        );
    } else {
        println!(
            "  {}",
            style("No CUDA detected — only the CPU option will be available").yellow()
        );
    }

    let mut labels: Vec<String> = Vec::new();
    let mut specs: Vec<String> = Vec::new();
    if cuda_count >= 1 {
        labels.push("Single GPU (cuda:0)".to_string());
        specs.push("cuda:0".to_string());
    }
    if cuda_count >= 2 {
        let all: Vec<String> = (0..cuda_count).map(|i| format!("cuda:{i}")).collect();
        labels.push(format!("All {cuda_count} GPUs ({})", all.join(",")));
        specs.push(all.join(","));
    }
    labels.push("CPU only".to_string());
    specs.push("cpu".to_string());
    labels.push("Custom (type a device spec)".to_string());
    specs.push(String::new());

    let default_idx = if let Some(spec) = default_spec {
        specs
            .iter()
            .position(|s| s == spec)
            .unwrap_or(labels.len() - 1)
    } else if cuda_count >= 2 {
        1
    } else if cuda_count == 1 {
        0
    } else {
        labels.iter().position(|l| l.starts_with("CPU")).unwrap_or(0)
    };

    let idx = Select::with_theme(theme)
        .items(&labels)
        .default(default_idx)
        .interact()
        .context("select devices")?;
    let spec = if specs[idx].is_empty() {

        Input::with_theme(theme)
            .with_prompt("  devices (comma-separated, e.g. cuda:0,cuda:2)")
            .default(default_spec.unwrap_or("cuda:0").to_string())
            .interact_text()
            .context("read device spec")?
    } else {
        specs[idx].clone()
    };
    let devices = Device::parse_list(&spec).context("parse devices")?;
    Ok(devices)
}

fn step_backpressure(
    theme: &ColorfulTheme,
    step: u8,
    default_max: usize,
    default_retry_ms: u64,
    device_count: usize,
) -> Result<(usize, u64)> {
    println!();
    println!("{}", step_header(step, 7, "Backpressure"));
    println!(
        "  {}",
        style(format!(
            "max-inflight is enforced PER device. With {device_count} device(s) the total cap is {} × {device_count} = {}.",
            default_max,
            default_max * device_count
        ))
        .dim()
    );
    let max_str: String = Input::with_theme(theme)
        .with_prompt("  max-inflight (per device)")
        .default(default_max.to_string())
        .interact_text()
        .context("read max-inflight")?;
    let max: usize = max_str.parse().context("max-inflight must be u64")?;
    let retry_str: String = Input::with_theme(theme)
        .with_prompt("  retry-after-ms")
        .default(default_retry_ms.to_string())
        .interact_text()
        .context("read retry-after-ms")?;
    let retry: u64 = retry_str.parse().context("retry-after-ms must be u64")?;
    Ok((max, retry))
}

fn step_metrics(
    theme: &ColorfulTheme,
    step: u8,
    default_bind: Option<&str>,
) -> Result<Option<String>> {
    println!();
    println!(
        "{}",
        step_header(step, 7, "Prometheus /metrics endpoint (optional)")
    );
    let want = Confirm::with_theme(theme)
        .with_prompt("  Expose /metrics for scraping?")
        .default(default_bind.is_some())
        .interact()
        .context("confirm metrics")?;
    if !want {
        return Ok(None);
    }
    let v: String = Input::with_theme(theme)
        .with_prompt("  metrics bind")
        .default(default_bind.unwrap_or("127.0.0.1:9090").to_string())
        .interact_text()
        .context("read metrics bind")?;
    Ok(Some(v))
}

fn step_tls(
    theme: &ColorfulTheme,
    step: u8,
    bind: &str,
    default_cert: Option<&Path>,
    default_key: Option<&Path>,
    default_client_ca: Option<&Path>,
) -> Result<(Option<PathBuf>, Option<PathBuf>, Option<PathBuf>)> {
    println!();
    println!("{}", step_header(step, 7, "TLS + mTLS (on by default)"));
    println!(
        "  {}",
        style(
            "Encrypts the gRPC channel. mTLS additionally requires every client \
             to present a certificate signed by a trusted CA.",
        )
        .dim()
    );
    let want = Confirm::with_theme(theme)
        .with_prompt("  Enable TLS?")
        .default(true)
        .interact()
        .context("confirm TLS")?;
    if !want {
        println!(
            "  {}",
            style("(TLS disabled — server will start in plaintext)").yellow()
        );
        return Ok((None, None, None));
    }
    let want_mtls = Confirm::with_theme(theme)
        .with_prompt("  Require mTLS (clients must present a cert)?")
        .default(true)
        .interact()
        .context("confirm mTLS")?;

    let have_files = default_cert.is_some() && default_key.is_some();
    let generate = if have_files {
        !Confirm::with_theme(theme)
            .with_prompt("  Use the existing cert/key files?")
            .default(true)
            .interact()
            .context("confirm use-existing TLS")?
    } else {
        Confirm::with_theme(theme)
            .with_prompt("  Generate a self-signed PKI now? (No = point at existing files)")
            .default(true)
            .interact()
            .context("confirm generate PKI")?
    };

    if generate {
        let dir_str: String = Input::with_theme(theme)
            .with_prompt("  output dir for generated certs")
            .default("clair3_connect_secrets".to_string())
            .interact_text()
            .context("read PKI dir")?;
        let dir = PathBuf::from(shellexpand_tilde(&dir_str));
        let sans = tls_sans_for_bind(bind);
        println!(
            "  {} generating CA + server + client certs (SANs: {})…",
            style("·").dim(),
            sans.join(", ")
        );
        let pki = crate::provisioning::generate_pki(&dir, &sans).context("generate PKI")?;
        let domain = sans.first().map(String::as_str).unwrap_or("localhost");
        println!("  {} wrote PKI under {}", style("✓").green(), dir.display());
        println!(
            "      server identity: {} + {}",
            pki.server_cert.display(),
            pki.server_key.display()
        );
        if want_mtls {
            println!("  {}", style("📁 Copy these files to the CLIENT machine:").bold().cyan());
            println!("        {}   (the server CA)", pki.ca_cert.display());
            println!("        {}   (this client's cert)", pki.client_cert.display());
            println!("        {}   (this client's key)", pki.client_key.display());
            println!(
                "  {} then run the client with: --tls-server-ca ca.crt --tls-client-cert client.crt --tls-client-key client.key --tls-domain {}",
                style("→").cyan(),
                domain
            );
        } else {
            println!("  {}", style("📁 Copy this file to the CLIENT machine:").bold().cyan());
            println!("        {}   (the server CA)", pki.ca_cert.display());
            println!(
                "  {} then run the client with: --tls-server-ca ca.crt --tls-domain {}",
                style("→").cyan(),
                domain
            );
        }
        let client_ca = if want_mtls { Some(pki.ca_cert) } else { None };
        Ok((Some(pki.server_cert), Some(pki.server_key), client_ca))
    } else {
        let cert = required_path(theme, "  TLS cert PEM (server identity)", default_cert)?;
        let key = required_path(theme, "  TLS key PEM (server private key)", default_key)?;
        let client_ca = if want_mtls {
            Some(required_path(
                theme,
                "  client-CA PEM (signs allowed client certs)",
                default_client_ca,
            )?)
        } else {
            None
        };
        Ok((Some(cert), Some(key), client_ca))
    }
}

fn tls_sans_for_bind(bind: &str) -> Vec<String> {
    let host = bind.rsplit_once(':').map(|(h, _)| h).unwrap_or(bind).trim();
    let mut sans = vec!["localhost".to_string(), "127.0.0.1".to_string()];
    if !host.is_empty()
        && host != "0.0.0.0"
        && host != "localhost"
        && host != "127.0.0.1"
    {
        sans.push(host.to_string());
    }
    sans
}

fn step_aead(
    theme: &ColorfulTheme,
    step: u8,
    default: Option<&Path>,
) -> Result<Option<PathBuf>> {
    println!();
    println!(
        "{}",
        step_header(step, 7, "Payload AEAD — AES-256-GCM (on by default)")
    );
    println!(
        "  {}",
        style(
            "Encrypts request/response payloads end-to-end, defending them even \
             from a compromised TLS-terminating sidecar. The same key must be on \
             every client.",
        )
        .dim()
    );
    let want = Confirm::with_theme(theme)
        .with_prompt("  Enable AEAD?")
        .default(true)
        .interact()
        .context("confirm AEAD")?;
    if !want {
        return Ok(None);
    }

    let generate = if default.is_some() {
        !Confirm::with_theme(theme)
            .with_prompt("  Use the existing key file?")
            .default(true)
            .interact()
            .context("confirm use-existing AEAD")?
    } else {
        Confirm::with_theme(theme)
            .with_prompt("  Generate a new 32-byte key now? (No = point at an existing file)")
            .default(true)
            .interact()
            .context("confirm generate AEAD")?
    };

    if generate {
        let path_str: String = Input::with_theme(theme)
            .with_prompt("  output path for the generated key")
            .default("clair3_connect_secrets/aead.key".to_string())
            .interact_text()
            .context("read AEAD key path")?;
        let path = PathBuf::from(shellexpand_tilde(&path_str));
        crate::provisioning::generate_aead_key(&path).context("generate AEAD key")?;
        println!(
            "  {} wrote a 32-byte key to {} (chmod 600)",
            style("✓").green(),
            path.display()
        );
        println!(
            "  {} {}",
            style("📁 Copy this file to the CLIENT machine too:").bold().cyan(),
            path.display()
        );
        println!(
            "  {} then run the client with: --aead-key {}",
            style("→").cyan(),
            path.display()
        );
        return Ok(Some(path));
    }

    let p = required_path(theme, "  AEAD key file (32 raw bytes)", default)?;

    if let Ok(meta) = std::fs::metadata(&p) {
        if meta.len() != 32 {
            println!(
                "  {}: file is {} bytes, expected 32. Generate with:",
                style("warning").yellow(),
                meta.len()
            );
            println!("    head -c 32 /dev/urandom > {}", p.display());
        }
    }
    Ok(Some(p))
}

fn required_path(
    theme: &ColorfulTheme,
    prompt: &str,
    default: Option<&Path>,
) -> Result<PathBuf> {
    let default_str = default
        .and_then(|p| p.to_str())
        .map(|s| s.to_string())
        .unwrap_or_default();
    loop {
        let s: String = Input::with_theme(theme)
            .with_prompt(prompt)
            .default(default_str.clone())
            .allow_empty(default_str.is_empty())
            .interact_text()
            .context("read path")?;
        if s.is_empty() {
            println!("  {}: path is required", style("error").red());
            continue;
        }
        let p = PathBuf::from(shellexpand_tilde(&s));
        if !p.exists() {
            println!(
                "  {}: file does not exist — try again or Ctrl-C to abort",
                style("error").red()
            );
            continue;
        }
        return Ok(p);
    }
}

fn print_summary(term: &Term, args: &InteractiveServerArgs) -> Result<()> {
    let line = "─".repeat(56);
    println!();
    term.write_line(&format!("{}", style(&line).cyan()))?;
    term.write_line(&format!("{}", style("  Summary").bold().cyan()))?;
    term.write_line(&format!("{}", style(&line).cyan()))?;
    let kv = |k: &str, v: &str| {
        println!("  {:<14} {}", style(k).dim(), v);
    };
    kv("Bind:", &args.bind);
    kv("Registry:", &args.registry.display().to_string());
    let dev_str = args
        .devices
        .iter()
        .map(|d| format!("{d}"))
        .collect::<Vec<_>>()
        .join(", ");
    kv("Devices:", &dev_str);
    kv(
        "Per-dev cap:",
        &format!(
            "{} (total {} across {} devices)",
            args.max_inflight,
            args.max_inflight * args.devices.len(),
            args.devices.len()
        ),
    );
    kv("Retry-after:", &format!("{} ms", args.retry_after_ms));
    kv(
        "Metrics:",
        &args
            .metrics_bind
            .clone()
            .unwrap_or_else(|| "off".to_string()),
    );
    kv(
        "TLS:",
        &match (&args.tls_cert, &args.tls_client_ca) {
            (Some(c), Some(ca)) => format!(
                "mTLS (cert {}, client-CA {})",
                c.display(),
                ca.display()
            ),
            (Some(c), None) => format!("one-way (cert {})", c.display()),
            (None, _) => "off (plaintext)".to_string(),
        },
    );
    kv(
        "AEAD:",
        &args
            .aead_key_path
            .as_ref()
            .map(|p| format!("AES-256-GCM (key {})", p.display()))
            .unwrap_or_else(|| "off".to_string()),
    );
    println!();
    Ok(())
}

fn shellexpand_tilde(s: &str) -> String {
    if let Some(rest) = s.strip_prefix('~') {
        if let Some(home) = std::env::var_os("HOME") {
            let mut p = PathBuf::from(home);
            if !rest.is_empty() {
                p.push(rest.trim_start_matches('/'));
            }
            return p.display().to_string();
        }
    }
    s.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tilde_expansion() {
        std::env::set_var("HOME", "/h/x");
        assert_eq!(shellexpand_tilde("~/foo.bam"), "/h/x/foo.bam");
        assert_eq!(shellexpand_tilde("~"), "/h/x");
        assert_eq!(shellexpand_tilde("/abs"), "/abs");
    }
}
