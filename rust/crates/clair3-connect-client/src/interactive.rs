use anyhow::{bail, Context, Result};
use console::{style, Term};
use dialoguer::{theme::ColorfulTheme, Confirm, Input, Select};
use std::path::{Path, PathBuf};

#[allow(missing_docs)]
#[derive(Debug, Clone)]
pub struct InteractiveArgs {
    pub bam: PathBuf,
    pub fasta: PathBuf,
    pub region: String,

    pub full_pipeline: bool,
    pub server: String,
    pub pileup_model: String,
    pub fa_model: String,
    pub platform: String,
    pub output: PathBuf,
    pub threads: usize,
    pub bed: Option<PathBuf>,

    pub tls_server_ca: Option<PathBuf>,

    pub tls_client_cert: Option<PathBuf>,

    pub tls_client_key: Option<PathBuf>,

    pub tls_domain: String,

    pub aead_key_path: Option<PathBuf>,
}

#[derive(Debug, Clone, Default)]
struct SecuritySettings {
    tls_server_ca: Option<PathBuf>,
    tls_client_cert: Option<PathBuf>,
    tls_client_key: Option<PathBuf>,
    tls_domain: String,
    aead_key_path: Option<PathBuf>,
}

pub fn run(defaults: InteractiveDefaults) -> Result<InteractiveArgs> {
    let theme = ColorfulTheme::default();
    let term = Term::stdout();

    print_banner(&term)?;

    let bam = step_path(
        &theme,
        1,
        "Indexed BAM (.bam + .bai/.csi)",
        defaults.bam.as_deref(),
        Some(check_bam),
    )?;

    let fasta = step_path(
        &theme,
        2,
        "Indexed reference FASTA (.fa + .fai)",
        defaults.fasta.as_deref(),
        Some(check_fasta),
    )?;

    let region = step_region(&theme, 3, &fasta, defaults.region.as_deref())?;

    let full_pipeline = step_mode(&theme, 4, defaults.full_pipeline)?;

    let (server, sec, probed) = connect_loop(&theme, 5, &defaults)?;

    let platform = step_platform(&theme, 6, defaults.platform.as_deref())?;

    let (pileup_model, fa_model) = pick_models(
        &theme,
        7,
        probed.as_ref(),
        defaults.pileup_model.as_deref(),
        defaults.fa_model.as_deref(),
        &platform,
        full_pipeline,
    )?;

    let threads = step_threads(&theme, 8, defaults.threads)?;

    let output = step_output(&theme, 9, defaults.output.as_deref())?;

    let bed = step_bed(&theme, defaults.bed.as_deref())?;

    let args = InteractiveArgs {
        bam,
        fasta,
        region,
        full_pipeline,
        server,
        pileup_model,
        fa_model,
        platform,
        output,
        threads,
        bed,
        tls_server_ca: sec.tls_server_ca,
        tls_client_cert: sec.tls_client_cert,
        tls_client_key: sec.tls_client_key,
        tls_domain: sec.tls_domain,
        aead_key_path: sec.aead_key_path,
    };

    print_summary(&term, &args)?;
    if !Confirm::with_theme(&theme)
        .with_prompt("Proceed with these settings?")
        .default(true)
        .interact()?
    {
        bail!("setup cancelled by user");
    }

    Ok(args)
}

#[derive(Debug, Clone, Default)]
pub struct InteractiveDefaults {

    pub bam: Option<PathBuf>,

    pub fasta: Option<PathBuf>,

    pub region: Option<String>,

    pub full_pipeline: bool,

    pub server: Option<String>,

    pub pileup_model: Option<String>,

    pub fa_model: Option<String>,

    pub platform: Option<String>,

    pub output: Option<PathBuf>,

    pub threads: Option<usize>,

    pub bed: Option<PathBuf>,

    pub tls_server_ca: Option<PathBuf>,

    pub tls_client_cert: Option<PathBuf>,

    pub tls_client_key: Option<PathBuf>,

    pub tls_domain: Option<String>,

    pub aead_key_path: Option<PathBuf>,
}

fn print_banner(term: &Term) -> Result<()> {
    let line = "─".repeat(56);
    term.write_line("")?;
    term.write_line(&format!("{}", style(&line).cyan()))?;
    term.write_line(&format!(
        "{}",
        style("  Clair3-connect — interactive setup").bold().cyan()
    ))?;
    term.write_line(&format!(
        "{}",
        style("  Re-run with --non-interactive to skip prompts").dim()
    ))?;
    term.write_line(&format!("{}", style(&line).cyan()))?;
    term.write_line("")?;
    Ok(())
}

fn step_header(step: u8, prompt: &str) -> String {
    format!(
        "{} {}",
        style(format!("[Step {step}/9]")).bold().green(),
        style(prompt).bold()
    )
}

fn step_mode(theme: &ColorfulTheme, step: u8, default_full: bool) -> Result<bool> {
    println!();
    println!("{}", step_header(step, "Calling mode"));
    println!(
        "  {}",
        style(
            "FA-only (default): only the full-alignment model runs over the \
             region(s) — no pileup model, no phasing, lower memory."
        )
        .dim()
    );
    println!(
        "  {}",
        style(
            "Full pipeline: classic Clair3 flow (pileup → SelectCandidates → \
             full-alignment → merge); supports --threads and --enable-phasing."
        )
        .dim()
    );
    let items = ["FA-only (default)", "Full pipeline (pileup + FA + merge)"];
    let idx = Select::with_theme(theme)
        .items(&items)
        .default(if default_full { 1 } else { 0 })
        .interact()
        .context("select calling mode")?;
    Ok(idx == 1)
}

fn step_path(
    theme: &ColorfulTheme,
    step: u8,
    prompt: &str,
    default: Option<&Path>,
    validator: Option<fn(&Path) -> Result<()>>,
) -> Result<PathBuf> {
    println!();
    println!("{}", step_header(step, prompt));
    let default_str: String = default
        .and_then(|p| p.to_str())
        .map(|s| s.to_string())
        .unwrap_or_default();
    loop {
        let input: String = Input::with_theme(theme)
            .with_prompt("  path")
            .default(default_str.clone())
            .interact_text()
            .context("read path")?;
        let path = PathBuf::from(shellexpand_tilde(&input));
        if !path.exists() {
            println!(
                "  {}: file does not exist — try again or Ctrl-C to abort",
                style("error").red()
            );
            continue;
        }
        if let Some(check) = validator {
            if let Err(e) = check(&path) {
                println!("  {}: {e}", style("error").red());
                continue;
            }
        }
        return Ok(path);
    }
}

fn check_bam(p: &Path) -> Result<()> {
    if p.extension().and_then(|s| s.to_str()).is_none_or(|s| s != "bam") {
        bail!("expected a .bam file");
    }
    let bai = p.with_extension("bam.bai");
    let csi = p.with_extension("bam.csi");
    if !bai.exists() && !csi.exists() {
        bail!("missing index — expected {}.bai or {}.csi alongside", p.display(), p.display());
    }
    Ok(())
}

fn check_fasta(p: &Path) -> Result<()> {
    let ok_ext = matches!(
        p.extension().and_then(|s| s.to_str()),
        Some("fa") | Some("fasta") | Some("fna")
    );
    if !ok_ext {
        bail!("expected .fa / .fasta / .fna");
    }
    let fai = PathBuf::from(format!("{}.fai", p.display()));
    if !fai.exists() {
        bail!("missing .fai index — run `samtools faidx {}`", p.display());
    }
    Ok(())
}

fn step_text(
    theme: &ColorfulTheme,
    step: u8,
    prompt: &str,
    default: &str,
) -> Result<String> {
    println!();
    println!("{}", step_header(step, prompt));
    let value: String = Input::with_theme(theme)
        .with_prompt("  value")
        .default(default.to_string())
        .interact_text()
        .with_context(|| format!("read {prompt}"))?;
    Ok(value)
}

fn step_region(
    theme: &ColorfulTheme,
    step: u8,
    fasta: &Path,
    default: Option<&str>,
) -> Result<String> {
    println!();
    println!("{}", step_header(step, "Region to call variants in"));
    println!(
        "  {}",
        style("Format: chrom:start-end (1-based inclusive). Examples:").dim()
    );
    println!(
        "    {}    {}",
        style("chr20:100000-300000").cyan(),
        style("(Phase-0 demo region)").dim()
    );
    println!(
        "    {}      {}",
        style("chr20:1-64444167").cyan(),
        style("(whole chr20)").dim()
    );

    let fai_default = first_contig_full_length(fasta).ok();
    let default_str = default
        .map(String::from)
        .or(fai_default)
        .unwrap_or_else(|| "chr20:100000-300000".to_string());
    let value: String = Input::with_theme(theme)
        .with_prompt("  region")
        .default(default_str)
        .interact_text()
        .context("read region")?;
    if !value.contains(':') || !value.contains('-') {
        bail!("region must be chrom:start-end, got {value:?}");
    }
    Ok(value)
}

fn first_contig_full_length(fasta: &Path) -> Result<String> {
    let fai = PathBuf::from(format!("{}.fai", fasta.display()));
    let text = std::fs::read_to_string(&fai).with_context(|| format!("read {}", fai.display()))?;
    let line = text.lines().next().context("empty .fai")?;
    let mut cols = line.split('\t');
    let name = cols.next().context("no contig name in .fai")?;
    let length: u64 = cols
        .next()
        .context("no length in .fai")?
        .parse()
        .context("parse contig length")?;
    Ok(format!("{name}:1-{length}"))
}

fn step_platform(
    theme: &ColorfulTheme,
    step: u8,
    default: Option<&str>,
) -> Result<String> {
    println!();
    println!("{}", step_header(step, "Sequencing platform"));
    let items = ["ont (Oxford Nanopore)", "hifi (PacBio HiFi)", "ilmn (Illumina)"];
    let default_idx = match default {
        Some("ont") => 0,
        Some("hifi") => 1,
        Some("ilmn") | Some("illumina") => 2,
        _ => 0,
    };
    let idx = Select::with_theme(theme)
        .items(&items)
        .default(default_idx)
        .interact()
        .context("select platform")?;
    Ok(match idx {
        0 => "ont".to_string(),
        1 => "hifi".to_string(),
        2 => "ilmn".to_string(),
        _ => unreachable!(),
    })
}

fn step_threads(
    theme: &ColorfulTheme,
    step: u8,
    default: Option<usize>,
) -> Result<usize> {
    println!();
    println!("{}", step_header(step, "Parallel chunks (threads)"));
    println!(
        "  {}",
        style("More threads = faster, more server pressure. 1 disables parallel driver.").dim()
    );
    let cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);
    let default = default.unwrap_or_else(|| std::cmp::min(cores, 8));
    let value: String = Input::with_theme(theme)
        .with_prompt("  threads")
        .default(default.to_string())
        .interact_text()
        .context("read threads")?;
    let n: usize = value.parse().context("threads must be a positive integer")?;
    if n == 0 {
        bail!("threads must be ≥ 1");
    }
    Ok(n)
}

fn step_output(
    theme: &ColorfulTheme,
    step: u8,
    default: Option<&Path>,
) -> Result<PathBuf> {
    println!();
    println!("{}", step_header(step, "Output directory"));
    let default_str = default
        .and_then(|p| p.to_str())
        .map(|s| s.to_string())
        .unwrap_or_else(|| "./clair3_connect_out".to_string());
    let value: String = Input::with_theme(theme)
        .with_prompt("  output dir")
        .default(default_str)
        .interact_text()
        .context("read output dir")?;
    Ok(PathBuf::from(shellexpand_tilde(&value)))
}

fn step_bed(theme: &ColorfulTheme, default: Option<&Path>) -> Result<Option<PathBuf>> {
    println!();
    println!(
        "{} {}",
        style("[Optional]").bold().yellow(),
        style("BED file to restrict variant calling to specific regions").bold()
    );
    let use_bed = Confirm::with_theme(theme)
        .with_prompt("  Use a BED inclusion file?")
        .default(default.is_some())
        .interact()
        .context("confirm BED")?;
    if !use_bed {
        return Ok(None);
    }
    let default_str = default
        .and_then(|p| p.to_str())
        .map(|s| s.to_string())
        .unwrap_or_default();
    loop {
        let input: String = Input::with_theme(theme)
            .with_prompt("  BED path")
            .default(default_str.clone())
            .interact_text()
            .context("read BED path")?;
        let path = PathBuf::from(shellexpand_tilde(&input));
        if !path.exists() {
            println!(
                "  {}: BED file does not exist — try again or skip",
                style("error").red()
            );
            continue;
        }
        return Ok(Some(path));
    }
}

fn run_async_in_existing_runtime<F, T>(fut: F) -> Result<T>
where
    F: std::future::Future<Output = Result<T>>,
{
    match tokio::runtime::Handle::try_current() {
        Ok(handle) => {

            Ok(tokio::task::block_in_place(move || handle.block_on(fut))?)
        }
        Err(_) => {

            let rt = tokio::runtime::Builder::new_current_thread()
                .enable_all()
                .build()?;
            rt.block_on(fut)
        }
    }
}

fn connect_loop(
    theme: &ColorfulTheme,
    step: u8,
    defaults: &InteractiveDefaults,
) -> Result<(String, SecuritySettings, Option<std::collections::HashMap<String, String>>)> {
    let mut server = defaults
        .server
        .clone()
        .unwrap_or_else(|| "http://127.0.0.1:50051".to_string());
    let mut sec = SecuritySettings {
        tls_server_ca: defaults.tls_server_ca.clone(),
        tls_client_cert: defaults.tls_client_cert.clone(),
        tls_client_key: defaults.tls_client_key.clone(),
        tls_domain: defaults
            .tls_domain
            .clone()
            .unwrap_or_else(|| "localhost".to_string()),
        aead_key_path: defaults.aead_key_path.clone(),
    };
    loop {
        println!();
        println!("{}", step_header(step, "Server endpoint + security"));
        server = prompt_text(theme, "gRPC server endpoint (http:// or https://)", &server)?;
        sec = step_security(theme, &sec)?;
        println!(
            "  {}",
            style(format!("Checking connection to {server} (Health RPC)…")).dim()
        );
        match probe_sync(&server, &sec) {
            Ok(map) => {
                println!(
                    "  {} connected — {} model(s) registered on the server",
                    style("✓").green(),
                    map.len()
                );
                return Ok((server, sec, Some(map)));
            }
            Err(e) => {
                print_conn_diagnosis(&sec, &e);
                let choice = Select::with_theme(theme)
                    .with_prompt("  What now?")
                    .items(&[
                        "Edit connection settings and retry",
                        "Continue anyway (type the model_id manually)",
                        "Quit",
                    ])
                    .default(0)
                    .interact()?;
                match choice {
                    0 => continue,
                    1 => return Ok((server, sec, None)),
                    _ => bail!("setup cancelled at the connection step"),
                }
            }
        }
    }
}

fn step_security(theme: &ColorfulTheme, cur: &SecuritySettings) -> Result<SecuritySettings> {
    let want_tls = Confirm::with_theme(theme)
        .with_prompt("  Is the server secured with TLS? (the server wizard defaults TLS ON)")
        .default(cur.tls_server_ca.is_some())
        .interact()?;
    if !want_tls {

        let aead = ask_optional_file(theme, "AEAD key file (32-byte; blank = none)", cur.aead_key_path.as_deref())?;
        return Ok(SecuritySettings {
            tls_server_ca: None,
            tls_client_cert: None,
            tls_client_key: None,
            tls_domain: cur.tls_domain.clone(),
            aead_key_path: aead,
        });
    }
    let ca = require_file(theme, "Server CA cert (ca.crt)", cur.tls_server_ca.as_deref())?;
    let mtls = Confirm::with_theme(theme)
        .with_prompt("  Does it require a client certificate (mTLS)? (server wizard defaults mTLS ON)")
        .default(cur.tls_client_cert.is_some())
        .interact()?;
    let (cert, key) = if mtls {
        (
            Some(require_file(theme, "Client cert (client.crt)", cur.tls_client_cert.as_deref())?),
            Some(require_file(theme, "Client key (client.key)", cur.tls_client_key.as_deref())?),
        )
    } else {
        (None, None)
    };
    let dflt_domain = if cur.tls_domain.is_empty() { "localhost".to_string() } else { cur.tls_domain.clone() };
    let domain = prompt_text(theme, "Server hostname on its cert (SNI)", &dflt_domain)?;
    let aead = ask_optional_file(theme, "AEAD key file (32-byte; blank = none)", cur.aead_key_path.as_deref())?;
    Ok(SecuritySettings {
        tls_server_ca: Some(ca),
        tls_client_cert: cert,
        tls_client_key: key,
        tls_domain: domain,
        aead_key_path: aead,
    })
}

fn prompt_text(theme: &ColorfulTheme, label: &str, initial: &str) -> Result<String> {
    Ok(Input::with_theme(theme)
        .with_prompt(format!("  {label}"))
        .with_initial_text(initial)
        .interact_text()?)
}

fn require_file(theme: &ColorfulTheme, label: &str, default: Option<&Path>) -> Result<PathBuf> {
    let init = default.map(|p| p.display().to_string()).unwrap_or_default();
    loop {
        let s: String = Input::with_theme(theme)
            .with_prompt(format!("  {label}"))
            .with_initial_text(&init)
            .interact_text()?;
        let p = PathBuf::from(shellexpand_tilde(s.trim()));
        if p.is_file() {
            return Ok(p);
        }
        println!("  {}: file not found — try again", style("error").red());
    }
}

fn ask_optional_file(theme: &ColorfulTheme, label: &str, default: Option<&Path>) -> Result<Option<PathBuf>> {
    let init = default.map(|p| p.display().to_string()).unwrap_or_default();
    loop {
        let s: String = Input::with_theme(theme)
            .with_prompt(format!("  {label}"))
            .with_initial_text(&init)
            .allow_empty(true)
            .interact_text()?;
        let s = s.trim();
        if s.is_empty() {
            return Ok(None);
        }
        let p = PathBuf::from(shellexpand_tilde(s));
        if !p.is_file() {
            println!("  {}: file not found — try again or leave blank", style("error").red());
            continue;
        }
        if let Ok(meta) = std::fs::metadata(&p) {
            if meta.len() != 32 {
                println!(
                    "  {}: an AEAD key must be exactly 32 bytes (this file is {})",
                    style("warning").yellow(),
                    meta.len()
                );
            }
        }
        return Ok(Some(p));
    }
}

fn probe_sync(
    server: &str,
    sec: &SecuritySettings,
) -> Result<std::collections::HashMap<String, String>> {
    run_async_in_existing_runtime(probe_server(server, sec))
}

async fn probe_server(
    server: &str,
    sec: &SecuritySettings,
) -> Result<std::collections::HashMap<String, String>> {
    use clair3_connect_transport::predict_client::{PredictClient, RetryPolicy};
    let policy = RetryPolicy {
        max_retries: 0,
        initial_delay: std::time::Duration::from_millis(100),
        max_delay: std::time::Duration::from_secs(1),
        per_attempt_timeout: Some(std::time::Duration::from_secs(5)),
    };
    let mut client = match sec.tls_server_ca.as_ref() {
        Some(ca) => {
            let tls = clair3_connect_transport::tls::ClientTlsOptions::from_paths(
                ca,
                sec.tls_client_cert.as_deref(),
                sec.tls_client_key.as_deref(),
                if sec.tls_domain.is_empty() { "localhost" } else { &sec.tls_domain },
            )
            .context("load TLS options")?;
            PredictClient::connect_tls_with_retry(server.to_string(), &tls, policy)
                .await
                .with_context(|| format!("TLS connect {server}"))?
        }
        None => PredictClient::connect_with_retry(server.to_string(), policy)
            .await
            .with_context(|| format!("connect {server}"))?,
    };
    let h = client.health().await.context("Health RPC")?;
    Ok(h.registered_models)
}

fn print_conn_diagnosis(sec: &SecuritySettings, err: &anyhow::Error) {
    let msg = format!("{err:#}");
    println!("  {} {}", style("connection failed:").red(), style(&msg).dim());
    let looks_h2 = msg.contains("h2 protocol") || msg.contains("http2") || msg.contains("frame");
    if sec.tls_server_ca.is_none() && (looks_h2 || msg.contains("transport")) {
        println!("  {}", style("→ The server appears to require TLS, but you connected in plaintext.").yellow());
        println!("    Answer YES to the TLS question and give the CA (ca.crt), the client");
        println!("    cert/key (client.crt/.key) if it uses mTLS, and the AEAD key if enabled.");
        println!("    The server's setup wizard prints exactly which files to copy here.");
    } else if sec.tls_server_ca.is_some() {
        let dom = if sec.tls_domain.is_empty() { "localhost" } else { &sec.tls_domain };
        println!("  {}", style("→ TLS is on but the handshake failed. Check:").yellow());
        println!("    • the CA / client cert+key match the server's PKI;");
        println!("    • the server hostname (SNI) matches its cert (this run: '{dom}');");
        println!("    • or the server is actually plaintext (answer NO to TLS, use http://).");
    } else {
        println!("  {}", style("→ Is the server running and is the endpoint correct?").yellow());
    }
}

#[allow(clippy::too_many_arguments)]
fn pick_models(
    theme: &ColorfulTheme,
    step: u8,
    probed: Option<&std::collections::HashMap<String, String>>,
    default_pileup: Option<&str>,
    default_fa: Option<&str>,
    platform: &str,
    full_pipeline: bool,
) -> Result<(String, String)> {
    println!();
    println!("{}", step_header(step, "Models"));
    let (pileup, fa) = if let Some(map) = probed {
        let mut ids: Vec<String> = map.keys().cloned().collect();
        ids.sort();
        if ids.is_empty() {
            println!("  {}: server has no registered models", style("error").red());
            bail!("server has no models");
        }

        let pileup_ids: Vec<String> = ids
            .iter()
            .filter(|s| s.contains("pileup"))
            .cloned()
            .collect();
        let fa_ids: Vec<String> = ids
            .iter()
            .filter(|s| s.contains("full_alignment") || s.starts_with("fa_"))
            .cloned()
            .collect();
        let f_choices: Vec<String> = if fa_ids.is_empty() { ids.clone() } else { fa_ids };

        let pileup = if full_pipeline {
            let p_choices: Vec<String> =
                if pileup_ids.is_empty() { ids.clone() } else { pileup_ids };
            println!();
            println!("  {}", style("Select pileup model:").bold());
            let p_default = p_choices
                .iter()
                .position(|s| Some(s.as_str()) == default_pileup)
                .or_else(|| p_choices.iter().position(|s| s.contains(platform)))
                .unwrap_or(0);
            let p_idx = Select::with_theme(theme)
                .items(&p_choices)
                .default(p_default)
                .interact()
                .context("select pileup model")?;
            p_choices[p_idx].clone()
        } else {
            pileup_ids
                .first()
                .cloned()
                .or_else(|| default_pileup.map(String::from))
                .unwrap_or_else(|| "pileup_ont_hac_v500".to_string())
        };

        println!("  {}", style("Select full-alignment model:").bold());
        let f_default = f_choices
            .iter()
            .position(|s| Some(s.as_str()) == default_fa)
            .or_else(|| f_choices.iter().position(|s| s.contains(platform)))
            .unwrap_or(0);
        let f_idx = Select::with_theme(theme)
            .items(&f_choices)
            .default(f_default)
            .interact()
            .context("select FA model")?;
        (pileup, f_choices[f_idx].clone())
    } else {
        let pileup = if full_pipeline {
            step_text(
                theme,
                step,
                "Pileup model_id",
                default_pileup.unwrap_or("pileup_ont_hac_v500"),
            )?
        } else {
            default_pileup.unwrap_or("pileup_ont_hac_v500").to_string()
        };
        let fa = step_text(
            theme,
            step,
            "Full-alignment model_id",
            default_fa.unwrap_or("full_alignment_ont_hac_v500"),
        )?;
        (pileup, fa)
    };
    Ok((pileup, fa))
}

fn print_summary(term: &Term, args: &InteractiveArgs) -> Result<()> {
    let line = "─".repeat(56);
    println!();
    term.write_line(&format!("{}", style(&line).cyan()))?;
    term.write_line(&format!("{}", style("  Summary").bold().cyan()))?;
    term.write_line(&format!("{}", style(&line).cyan()))?;
    let kv = |k: &str, v: &str| {
        println!("  {:<14} {}", style(k).dim(), v);
    };
    kv("BAM:", &args.bam.display().to_string());
    kv("FASTA:", &args.fasta.display().to_string());
    kv("Region:", &args.region);
    if let Some(bed) = args.bed.as_ref() {
        kv("BED:", &bed.display().to_string());
    }
    kv(
        "Mode:",
        if args.full_pipeline {
            "full pipeline (pileup + FA + merge)"
        } else {
            "FA-only (default)"
        },
    );
    kv("Server:", &args.server);
    let security = match (&args.tls_server_ca, args.tls_client_cert.is_some(), &args.aead_key_path) {
        (Some(_), true, Some(_)) => "mTLS + AEAD".to_string(),
        (Some(_), true, None) => "mTLS".to_string(),
        (Some(_), false, Some(_)) => "TLS + AEAD".to_string(),
        (Some(_), false, None) => "TLS".to_string(),
        (None, _, Some(_)) => "AEAD only (plaintext channel)".to_string(),
        (None, _, None) => "plaintext".to_string(),
    };
    kv("Security:", &security);
    kv("Platform:", &args.platform);

    if args.full_pipeline {
        kv("Pileup model:", &args.pileup_model);
    }
    kv("FA model:", &args.fa_model);
    if args.full_pipeline {
        kv("Threads:", &args.threads.to_string());
    }
    kv("Output:", &args.output.display().to_string());
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
    fn tilde_expansion_replaces_home_prefix() {
        std::env::set_var("HOME", "/home/test");
        assert_eq!(shellexpand_tilde("~/x.bam"), "/home/test/x.bam");
        assert_eq!(shellexpand_tilde("~"), "/home/test");
        assert_eq!(shellexpand_tilde("/abs/x.bam"), "/abs/x.bam");
        assert_eq!(shellexpand_tilde("rel/x.bam"), "rel/x.bam");
    }

    fn unique_tmp(label: &str, ext: &str) -> PathBuf {
        static COUNTER: std::sync::atomic::AtomicU64 =
            std::sync::atomic::AtomicU64::new(0);
        let n = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "clair3-connect-iv-{}-{}-{}.{ext}",
            std::process::id(),
            label,
            n
        ))
    }

    #[test]
    fn check_bam_requires_index() {
        let tmp = unique_tmp("bam", "bam");
        std::fs::write(&tmp, b"").unwrap();
        let res = check_bam(&tmp);
        assert!(res.is_err());
        let bai = tmp.with_extension("bam.bai");
        std::fs::write(&bai, b"").unwrap();
        check_bam(&tmp).unwrap();
        std::fs::remove_file(&tmp).ok();
        std::fs::remove_file(&bai).ok();
    }

    #[test]
    fn check_fasta_requires_fai() {
        let tmp = unique_tmp("fa-noidx", "fa");
        std::fs::write(&tmp, b"").unwrap();
        assert!(check_fasta(&tmp).is_err());
        let fai = PathBuf::from(format!("{}.fai", tmp.display()));
        std::fs::write(&fai, b"chr1\t100\t6\t60\t61\n").unwrap();
        check_fasta(&tmp).unwrap();
        std::fs::remove_file(&tmp).ok();
        std::fs::remove_file(&fai).ok();
    }

    #[test]
    fn first_contig_full_length_reads_fai() {
        let tmp = unique_tmp("fa-len", "fa");
        std::fs::write(&tmp, b"").unwrap();
        let fai = PathBuf::from(format!("{}.fai", tmp.display()));
        std::fs::write(
            &fai,
            b"chr20\t64444167\t6\t60\t61\nchr21\t46709983\t...\n",
        )
        .unwrap();
        let r = first_contig_full_length(&tmp).unwrap();
        assert_eq!(r, "chr20:1-64444167");
        std::fs::remove_file(&tmp).ok();
        std::fs::remove_file(&fai).ok();
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn run_async_in_existing_runtime_does_not_panic() {

        let r = run_async_in_existing_runtime(async {
            Ok(("pileup".to_string(), "fa".to_string()))
        });
        let (p, f) = r.expect("must not panic");
        assert_eq!(p, "pileup");
        assert_eq!(f, "fa");
    }

    #[test]
    fn run_async_in_existing_runtime_works_outside_any_runtime() {
        let r = run_async_in_existing_runtime(async {
            Ok(("p".to_string(), "f".to_string()))
        });
        let (p, f) = r.expect("must not panic");
        assert_eq!(p, "p");
        assert_eq!(f, "f");
    }

    #[test]
    fn check_fasta_accepts_common_extensions() {
        for ext in ["fa", "fasta", "fna"] {
            let tmp = unique_tmp("fa-ext", ext);
            std::fs::write(&tmp, b"").unwrap();
            let fai = PathBuf::from(format!("{}.fai", tmp.display()));
            std::fs::write(&fai, b"chr1\t100\n").unwrap();
            check_fasta(&tmp).unwrap_or_else(|e| panic!("ext {ext}: {e}"));
            std::fs::remove_file(&tmp).ok();
            std::fs::remove_file(&fai).ok();
        }
    }
}
