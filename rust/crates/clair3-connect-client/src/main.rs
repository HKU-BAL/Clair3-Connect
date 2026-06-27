use anyhow::{Context, Result};
use clair3_connect_client::chunk_planner::ChunkRegion;
use clair3_connect_client::fa_only_pipeline::{run_fa_only_pipeline, FaOnlyInputs};
use clair3_connect_client::fa_pipeline::{run_fa_pipeline, FaPipelineInputs};
use clair3_connect_client::merge_vcf::{merge_vcf, MergeConfig};
use clair3_connect_client::parallel_driver::{run_parallel_pileup, ParallelPileupConfig};
use clair3_connect_client::pileup_pipeline::{
    expand_user_region, run_pileup_pipeline, PileupPipelineInputs,
};
use clair3_connect_client::region_filter::{filter_vcf_by_region, RegionFilter};
use clair3_connect_client::select_candidates::{select_candidates, CandidateFilterMode};
use clair3_connect_client::select_het_snp::{filter_het_snp_phasing, HetSnpFilterMode};
use clair3_connect_client::select_qual::{select_phase_qual, select_qual};
use clair3_connect_core::full_alignment::phase_block_index::load_phased_vcf;
use clair3_connect_core::hts::region::Region;
use clair3_connect_longphase::runner::{run_longphase_phase, LongPhaseInputs, LongPhasePlatform};
use clair3_connect_core::pileup::kernel::PileupParams;
use clair3_connect_core::platform::{Platform, PlatformParams};
use clair3_connect_core::vcf::output_with::OutputWithConfig;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Debug)]
struct Args {

    bam: Option<PathBuf>,
    fasta: Option<PathBuf>,
    region: Option<String>,

    bed: Option<PathBuf>,

    regions: Vec<String>,

    aead_key_path: Option<PathBuf>,

    tls_server_ca: Option<PathBuf>,

    tls_client_cert: Option<PathBuf>,

    tls_client_key: Option<PathBuf>,

    tls_domain: String,
    pileup_model: String,
    fa_model: String,
    server: String,
    output: PathBuf,
    sample_name: String,
    clair3_version: String,
    platform: String,
    var_pct_full: f64,
    ref_pct_full: f64,
    threads: usize,

    threads_set: bool,

    max_concurrent: Option<usize>,
    chunk_size: u64,

    low_mem: bool,

    fa_window_bp: u64,

    full_pipeline: bool,

    enable_phasing: bool,

    longphase: PathBuf,

    var_pct_phasing: f64,

    force_interactive: bool,

    force_non_interactive: bool,
}

fn parse_args() -> Args {
    let mut bam: Option<PathBuf> = None;
    let mut fasta: Option<PathBuf> = None;

    let mut region_opt: Option<String> = None;
    let mut bed: Option<PathBuf> = None;
    let mut regions: Vec<String> = Vec::new();
    let mut aead_key_path: Option<PathBuf> = None;
    let mut tls_server_ca: Option<PathBuf> = None;
    let mut tls_client_cert: Option<PathBuf> = None;
    let mut tls_client_key: Option<PathBuf> = None;
    let mut tls_domain: String = "localhost".to_string();
    let mut pileup_model = "pileup_ont_hac_v500".to_string();
    let mut fa_model = "full_alignment_ont_hac_v500".to_string();
    let mut server = "http://127.0.0.1:50051".to_string();
    let mut output: PathBuf = PathBuf::from("clair3_connect_out");
    let mut sample_name = "SAMPLE".to_string();
    let mut clair3_version = "2.0.1".to_string();
    let mut platform = "ont".to_string();
    let mut var_pct_full = 0.3_f64;
    let mut ref_pct_full = 0.3_f64;
    let mut threads = 1_usize;
    let mut threads_set = false;
    let mut max_concurrent: Option<usize> = None;
    let mut chunk_size = 5_000_000_u64;
    let mut low_mem = false;
    let mut fa_window_bp = 0u64;

    let mut full_pipeline = false;
    let mut enable_phasing = false;
    let mut longphase = PathBuf::from("longphase");
    let mut var_pct_phasing = 0.7_f64;
    let mut force_interactive = false;
    let mut force_non_interactive = false;

    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        macro_rules! take {
            () => {
                args.next().unwrap_or_else(|| panic!("flag {} needs a value", a))
            };
        }
        match a.as_str() {
            "--bam" => bam = Some(PathBuf::from(take!())),
            "--ref" => fasta = Some(PathBuf::from(take!())),
            "--region" => region_opt = Some(take!()),
            "--interactive" | "-i" => force_interactive = true,
            "--non-interactive" => force_non_interactive = true,
            "--bed" => bed = Some(PathBuf::from(take!())),
            "--aead-key" => aead_key_path = Some(PathBuf::from(take!())),
            "--tls-server-ca" => tls_server_ca = Some(PathBuf::from(take!())),
            "--tls-client-cert" => tls_client_cert = Some(PathBuf::from(take!())),
            "--tls-client-key" => tls_client_key = Some(PathBuf::from(take!())),
            "--tls-domain" => tls_domain = take!(),
            "--regions" => {

                let s = take!();
                for piece in s.split(',') {
                    let t = piece.trim();
                    if !t.is_empty() {
                        regions.push(t.to_string());
                    }
                }
            }
            "--pileup-model" => pileup_model = take!(),
            "--fa-model" => fa_model = take!(),
            "--server" => server = take!(),
            "--output" => output = PathBuf::from(take!()),
            "--sample-name" => sample_name = take!(),
            "--clair3-version" => clair3_version = take!(),
            "--platform" => platform = take!(),
            "--var-pct-full" => var_pct_full = take!().parse().expect("--var-pct-full"),
            "--ref-pct-full" => ref_pct_full = take!().parse().expect("--ref-pct-full"),
            "--threads" => {
                threads = take!().parse().expect("--threads");
                threads_set = true;
            }
            "--max-concurrent" => {
                max_concurrent = Some(take!().parse().expect("--max-concurrent"))
            }
            "--chunk-size" => chunk_size = take!().parse().expect("--chunk-size"),
            "--low-mem" => low_mem = true,
            "--fa-window-bp" => fa_window_bp = take!().parse().expect("--fa-window-bp"),

            "--fa-only" => full_pipeline = false,
            "--full-pipeline" => full_pipeline = true,
            "--enable-phasing" => enable_phasing = true,
            "--longphase" => longphase = PathBuf::from(take!()),
            "--var-pct-phasing" => var_pct_phasing = take!().parse().expect("--var-pct-phasing"),
            "-h" | "--help" => {
                eprintln!(
                    "clair3-connect --bam BAM --ref REF.fa --region chr:start-end \\
        --fa-model ID --server http://host:port \\
        --output DIR [--sample-name S] [--clair3-version V] \\
        [--platform ont|hifi|ilmn] \\
        [--bed PATH] [--regions chr:start-end[,chr:start-end...]] [--fa-window-bp BP] \\
        [--aead-key PATH (32-byte binary)] \\
        [--tls-server-ca PEM [--tls-client-cert PEM --tls-client-key PEM] [--tls-domain HOST]] \\
        [--interactive | --non-interactive]

Default mode is FA-only: ONLY the full-alignment model is run over the
regions of interest (--region / --bed / --regions); no pileup model, no
phasing. Final calls land in merge_output.vcf. (--fa-only re-asserts this
default and is accepted for backward compatibility.)

FA-only builds windows in parallel by default (auto threads = min(cores, 8),
500 kb windows) — no tuning needed. Override with --threads N and
--fa-window-bp BP (peak RSS ~ threads x per-window; --threads 1 = serial).

Pass --full-pipeline for the classic Clair3 flow (pileup -> SelectQual ->
SelectCandidates -> full-alignment -> merge), which additionally honours:
        --pileup-model ID [--var-pct-full 0.3] [--ref-pct-full 0.3] \\
        [--threads 1] [--max-concurrent N] [--chunk-size 5000000] [--low-mem] \\
        [--enable-phasing [--longphase PATH] [--var-pct-phasing 0.7]]

Tip: run with no required flags (or with --interactive) to start a
guided setup that walks you through every input."
                );
                std::process::exit(0);
            }
            other => {
                eprintln!("unknown flag: {other}");
                std::process::exit(1);
            }
        }
    }
    if tls_client_cert.is_some() != tls_client_key.is_some() {
        eprintln!("--tls-client-cert and --tls-client-key must be specified together");
        std::process::exit(1);
    }
    Args {
        bam,
        fasta,
        region: region_opt,
        bed,
        regions,
        aead_key_path,
        tls_server_ca,
        tls_client_cert,
        tls_client_key,
        tls_domain,
        pileup_model,
        fa_model,
        server,
        output,
        sample_name,
        clair3_version,
        platform,
        var_pct_full,
        ref_pct_full,
        threads,
        threads_set,
        max_concurrent,
        chunk_size,
        low_mem,
        fa_window_bp,
        full_pipeline,
        enable_phasing,
        longphase,
        var_pct_phasing,
        force_interactive,
        force_non_interactive,
    }
}

fn should_run_interactive(args: &Args) -> bool {
    if args.force_non_interactive {
        return false;
    }
    if args.force_interactive {
        return true;
    }
    let missing_required =
        args.bam.is_none() || args.fasta.is_none() || args.region.is_none();
    let on_tty = console::user_attended();
    missing_required && on_tty
}

fn merge_interactive(args: &mut Args, iv: clair3_connect_client::interactive::InteractiveArgs) {
    if args.bam.is_none() {
        args.bam = Some(iv.bam);
    }
    if args.fasta.is_none() {
        args.fasta = Some(iv.fasta);
    }
    if args.region.is_none() {
        args.region = Some(iv.region);
    }
    if args.bed.is_none() {
        args.bed = iv.bed;
    }

    args.full_pipeline = iv.full_pipeline;

    if args.server == "http://127.0.0.1:50051" {
        args.server = iv.server;
    }
    if args.pileup_model == "pileup_ont_hac_v500" {
        args.pileup_model = iv.pileup_model;
    }
    if args.fa_model == "full_alignment_ont_hac_v500" {
        args.fa_model = iv.fa_model;
    }
    if args.platform == "ont" {
        args.platform = iv.platform;
    }
    if args.threads == 1 {
        args.threads = iv.threads;
    }
    if args.output == Path::new("clair3_connect_out") {
        args.output = iv.output;
    }

    if args.tls_server_ca.is_none() {
        args.tls_server_ca = iv.tls_server_ca;
    }
    if args.tls_client_cert.is_none() {
        args.tls_client_cert = iv.tls_client_cert;
    }
    if args.tls_client_key.is_none() {
        args.tls_client_key = iv.tls_client_key;
    }
    if args.tls_domain == "localhost" && !iv.tls_domain.is_empty() {
        args.tls_domain = iv.tls_domain;
    }
    if args.aead_key_path.is_none() {
        args.aead_key_path = iv.aead_key_path;
    }
}

fn build_region_filter(args: &Args) -> Result<RegionFilter> {
    let mut filter = RegionFilter::empty();
    if let Some(bed) = args.bed.as_ref() {
        let f = RegionFilter::from_bed_path(bed)
            .with_context(|| format!("parse BED {}", bed.display()))?;
        filter = filter.union(f);
    }
    if !args.regions.is_empty() {

        let lens: HashMap<String, u64> = HashMap::new();
        let specs: Vec<&str> = args.regions.iter().map(|s| s.as_str()).collect();
        let f = RegionFilter::from_region_strings(&specs, &lens)
            .context("parse --regions")?;
        filter = filter.union(f);
    }
    Ok(filter)
}

fn pparams_str(p: Platform) -> &'static str {
    p.as_str()
}

fn parse_region_arg(region: &str) -> (String, u64, u64) {
    let (contig, rng) = region.split_once(':').expect("region needs ':'");
    let (start_s, end_s) = rng.split_once('-').expect("region needs '-'");
    let start: u64 = start_s.parse().expect("region start");
    let end: u64 = end_s.parse().expect("region end");
    (contig.to_string(), start, end)
}

#[tokio::main]
async fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();

    let mut args = parse_args();

    if should_run_interactive(&args) {
        let defaults = clair3_connect_client::interactive::InteractiveDefaults {
            bam: args.bam.clone(),
            fasta: args.fasta.clone(),
            region: args.region.clone(),
            full_pipeline: args.full_pipeline,
            server: Some(args.server.clone()),
            pileup_model: Some(args.pileup_model.clone()),
            fa_model: Some(args.fa_model.clone()),
            platform: Some(args.platform.clone()),
            output: Some(args.output.clone()),
            threads: Some(args.threads),
            bed: args.bed.clone(),
            tls_server_ca: args.tls_server_ca.clone(),
            tls_client_cert: args.tls_client_cert.clone(),
            tls_client_key: args.tls_client_key.clone(),
            tls_domain: Some(args.tls_domain.clone()),
            aead_key_path: args.aead_key_path.clone(),
        };
        let iv = clair3_connect_client::interactive::run(defaults)
            .context("interactive setup")?;
        merge_interactive(&mut args, iv);
    }

    let bam = args
        .bam
        .clone()
        .context("--bam is required (or run interactively without --non-interactive)")?;
    let fasta = args
        .fasta
        .clone()
        .context("--ref is required (or run interactively without --non-interactive)")?;
    let region = args
        .region
        .clone()
        .unwrap_or_else(|| "chr20:100000-300000".to_string());

    if args.chunk_size < clair3_connect_client::chunk_planner::MIN_CHUNK_SIZE {
        anyhow::bail!(
            "--chunk-size {} below MIN_CHUNK_SIZE {}",
            args.chunk_size,
            clair3_connect_client::chunk_planner::MIN_CHUNK_SIZE
        );
    }
    std::fs::create_dir_all(&args.output).context("mkdir output")?;

    let tls: Option<clair3_connect_transport::tls::ClientTlsOptions> = match args.tls_server_ca.as_ref() {
        Some(ca_path) => {
            let opts = clair3_connect_transport::tls::ClientTlsOptions::from_paths(
                ca_path,
                args.tls_client_cert.as_deref(),
                args.tls_client_key.as_deref(),
                &args.tls_domain,
            )
            .context("load TLS options")?;
            tracing::info!(
                ca = %ca_path.display(),
                mtls = opts.client_identity_pem.is_some(),
                domain = %args.tls_domain,
                "TLS enabled on client"
            );
            Some(opts)
        }
        None => None,
    };

    let aead_key: Option<clair3_connect_transport::aead::AeadKey> = match args.aead_key_path.as_ref() {
        Some(p) => {
            let bytes = std::fs::read(p)
                .with_context(|| format!("read --aead-key {}", p.display()))?;
            if bytes.len() != clair3_connect_transport::aead::AEAD_KEY_BYTES {
                anyhow::bail!(
                    "--aead-key file {} must be {} bytes (got {})",
                    p.display(),
                    clair3_connect_transport::aead::AEAD_KEY_BYTES,
                    bytes.len()
                );
            }
            let k = clair3_connect_transport::aead::AeadKey::from_bytes(&bytes)
                .context("parse AEAD key")?;
            tracing::info!(path = %p.display(), "AEAD enabled (AES-256-GCM)");
            Some(k)
        }
        None => None,
    };

    let region_filter = build_region_filter(&args).context("build region filter")?;
    if !region_filter.is_unrestricted() {
        let total_bp: u64 = region_filter
            .contigs()
            .map(|c| region_filter.covered_bp(c))
            .sum();
        tracing::info!(
            "region filter active: {} contigs, {} bp total",
            region_filter.contigs().count(),
            total_bp
        );
    }

    let (contig, ctg_start, ctg_end) = parse_region_arg(&region);

    let platform = Platform::from_str_or_ont(args.platform.as_str());
    let pparams = PlatformParams::for_platform(platform);

    let max_depth: i64 = pparams.max_depth as i64;
    let matrix_depth: usize = pparams.matrix_depth;
    tracing::info!(
        platform = pparams_str(platform),
        max_depth, matrix_depth,
        min_af = pparams.min_af,
        "platform parameters resolved"
    );

    if !args.full_pipeline {
        let rois: Vec<Region> = if !region_filter.is_unrestricted() {
            let mut v = Vec::new();
            for c in region_filter.contigs() {
                if let Some(ivs) = region_filter.intervals_for(c) {
                    for iv in ivs {
                        v.push(Region {
                            contig: c.to_string(),
                            start: iv.start,
                            end: iv.end,
                        });
                    }
                }
            }
            v
        } else {

            vec![Region {
                contig: contig.clone(),
                start: ctg_start.saturating_sub(1),
                end: ctg_end,
            }]
        };
        let total_roi_bp: u64 = rois.iter().map(|r| r.end.saturating_sub(r.start)).sum();
        tracing::info!(
            "[fa-only] {} ROIs, {} bp total — no pileup model, no phasing",
            rois.len(),
            total_roi_bp
        );

        let fa_only_inputs = FaOnlyInputs {
            bam_path: &bam,
            fasta_path: &fasta,
            rois,
            server_endpoint: args.server.clone(),
            fa_model_id: args.fa_model.clone(),
            sample_name: args.sample_name.clone(),
            clair3_version: args.clair3_version.clone(),
            matrix_depth,
            max_indel_length: 50,
            min_mq: 5,
            pileup_params: PileupParams::default(),
            output_config: OutputWithConfig {
                is_show_reference: false,
                pileup: false,
                add_indel_length: true,

                quality_score_for_pass: Some(2.0),
                ..OutputWithConfig::default()
            },
            batch_size: 32,
            fa_window_bp: args.fa_window_bp,

            fa_threads: if args.threads_set {
                args.threads.max(1)
            } else {
                std::thread::available_parallelism()
                    .map(|n| n.get().min(8))
                    .unwrap_or(4)
            },
            aead_key: aead_key.clone(),
            tls: tls.clone(),
        };

        let mut fa_buf: Vec<u8> = Vec::new();
        run_fa_only_pipeline(&fa_only_inputs, &mut fa_buf)
            .await
            .context("fa-only pipeline")?;

        let fa_vcf_path = args.output.join("full_alignment.vcf");
        let merge_vcf_path = args.output.join("merge_output.vcf");
        for p in [&fa_vcf_path, &merge_vcf_path] {
            let mut f = BufWriter::new(File::create(p)?);
            f.write_all(&fa_buf)?;
            f.flush()?;
        }
        let n_variants = String::from_utf8_lossy(&fa_buf)
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
            .count();
        tracing::info!("[fa-only] wrote {}", merge_vcf_path.display());
        println!(
            "Wrote {} ({} variants, FA-only).",
            merge_vcf_path.display(),
            n_variants
        );
        return Ok(());
    }

    tracing::info!("[1/4] pileup pipeline (threads={})", args.threads);
    let pileup_vcf_path = args.output.join("pileup.vcf");
    let pileup_params = PileupParams::default();
    let pileup_output_config = OutputWithConfig {
        is_show_reference: true,
        ..OutputWithConfig::default()
    };
    let mut pileup_buf: Vec<u8> = Vec::new();
    if args.threads > 1 && ctg_end - ctg_start > args.chunk_size {
        let mut chunks: Vec<ChunkRegion> = Vec::new();
        let mut cs = ctg_start;
        let mut idx = 0_u32;
        let span = ctg_end - ctg_start;
        let chunk_count = span.div_ceil(args.chunk_size) as u32;
        while cs <= ctg_end {
            let ce = (cs + args.chunk_size - 1).min(ctg_end);
            chunks.push(ChunkRegion {
                contig: contig.clone(),
                start_1based: cs,
                end_1based: ce,
                chunk_idx: idx,
                chunk_count,
            });
            cs = ce + 1;
            idx += 1;
        }

        if !region_filter.is_unrestricted() {
            let before = chunks.len();
            chunks.retain(|c| {
                !region_filter
                    .intersect_chunk_1based(&c.contig, c.start_1based, c.end_1based)
                    .is_empty()
            });
            let skipped = before - chunks.len();
            if skipped > 0 {
                tracing::info!(
                    "region filter: skipped {} of {} pileup chunks (no BED overlap)",
                    skipped,
                    before
                );
            }
        }

        let resume_dir = args.output.join("resume_state");

        let eff_max_concurrent = args.max_concurrent.unwrap_or_else(|| {
            if args.low_mem { args.threads.min(4) } else { args.threads }
        });
        if args.low_mem && args.max_concurrent.is_none() && eff_max_concurrent < args.threads {
            tracing::info!(
                "low-mem: capping pileup max_concurrent to {} (of {} threads) to bound peak RSS",
                eff_max_concurrent, args.threads
            );
        }
        let par_config = ParallelPileupConfig {
            bam_path: bam.clone(),
            fasta_path: fasta.clone(),
            server_endpoint: args.server.clone(),
            model_id: args.pileup_model.clone(),
            sample_name: args.sample_name.clone(),
            clair3_version: args.clair3_version.clone(),
            max_depth,
            pileup_params: pileup_params.clone(),
            output_config: pileup_output_config,
            batch_size: 32,
            max_concurrent: eff_max_concurrent,
            manifest_dir: Some(resume_dir),
            verify_resume: true,
            aead_key: aead_key.clone(),
            tls: tls.clone(),
        };
        let stats = run_parallel_pileup(&chunks, &par_config, &mut pileup_buf)
            .await
            .context("parallel pileup")?;
        tracing::info!("parallel pileup: {} chunks, {} rows", stats.chunks_completed, stats.rows_emitted);
    } else {
        let pileup_inputs = PileupPipelineInputs {
            bam_path: &bam,
            fasta_path: &fasta,
            region: expand_user_region(&contig, ctg_start, ctg_end),
            server_endpoint: args.server.clone(),
            model_id: args.pileup_model.clone(),
            sample_name: args.sample_name.clone(),
            clair3_version: args.clair3_version.clone(),
            max_depth,
            pileup_params: pileup_params.clone(),
            output_config: pileup_output_config,
            batch_size: 32,
            aead_key: aead_key.clone(),
            tls: tls.clone(),
        };
        run_pileup_pipeline(&pileup_inputs, &mut pileup_buf)
            .await
            .context("pileup pipeline")?;
    }
    let pileup_text = String::from_utf8(pileup_buf).context("pileup utf8")?;
    let mut f = BufWriter::new(File::create(&pileup_vcf_path)?);
    f.write_all(pileup_text.as_bytes())?;
    f.flush()?;
    tracing::info!("[1/4] wrote {}", pileup_vcf_path.display());

    tracing::info!("[2/4] select_qual");
    let qual_sel =
        select_qual(pileup_text.as_bytes(), args.var_pct_full, args.ref_pct_full)
            .context("select_qual")?;
    tracing::info!(
        "var_cutoff={:.2}, ref_cutoff={:.2} ({} variants, {} refs)",
        qual_sel.var_qual_cut_off,
        qual_sel.ref_qual_cut_off,
        qual_sel.var_count,
        qual_sel.ref_count
    );

    tracing::info!("[3/4] select_candidates + FA pipeline");
    let cand_sel = select_candidates(
        pileup_text.as_bytes(),
        Some(&contig),
        CandidateFilterMode::QualCutoffs {
            var_cutoff: qual_sel.var_qual_cut_off,
            ref_cutoff: qual_sel.ref_qual_cut_off,
        },
    )
    .context("select_candidates")?;
    tracing::info!(
        "fa candidates: {} ({} variants, {} refs)",
        cand_sel.positions.len(),
        cand_sel.variants_kept,
        cand_sel.refs_kept
    );
    let candidates_0b: Vec<u64> = {
        let mut v: Vec<u64> = cand_sel.positions.iter().map(|p| p - 1).collect();
        v.sort();
        v.dedup();

        if !region_filter.is_unrestricted() {
            let before = v.len();
            v.retain(|p| region_filter.contains_position(&contig, *p));
            tracing::info!(
                "region filter: {} -> {} candidates ({} dropped)",
                before,
                v.len(),
                before - v.len()
            );
        }
        v
    };

    let phased_variants = if args.enable_phasing {
        let lp_platform = match args.platform.as_str() {
            "ont" => Some(LongPhasePlatform::Ont),
            "hifi" => Some(LongPhasePlatform::PacBio),
            other => {
                tracing::warn!(
                    "--enable-phasing: LongPhase is long-read only; skipping phasing for '{other}'"
                );
                None
            }
        };
        match lp_platform {
            Some(plat) => {
                tracing::info!("[phasing] LongPhase haplotagging enabled (platform={plat:?})");
                let phase_sel = select_phase_qual(pileup_text.as_bytes(), args.var_pct_phasing)
                    .context("select_phase_qual")?;
                let hetsnp_path = args.output.join("phasing_input_hetsnp.vcf");
                {
                    let mut w = BufWriter::new(File::create(&hetsnp_path)?);
                    filter_het_snp_phasing(
                        pileup_text.as_bytes(),
                        &mut w,
                        Some(contig.as_str()),
                        HetSnpFilterMode::QualCutoff(phase_sel.qual_cut_off),
                    )
                    .context("filter_het_snp_phasing")?;
                    w.flush()?;
                }
                let lp_prefix = args.output.join("phased");
                let lp_out = run_longphase_phase(&LongPhaseInputs {
                    binary: &args.longphase,
                    input_vcf: &hetsnp_path,
                    bam_path: &bam,
                    fasta_path: &fasta,
                    output_prefix: &lp_prefix,
                    threads: args.threads as u32,
                    platform: plat,
                })
                .context("run_longphase_phase")?;
                let variants = load_phased_vcf(&lp_out.phased_vcf, Some(&contig))
                    .context("load_phased_vcf")?;
                tracing::info!("[phasing] {} phased het SNPs from LongPhase", variants.len());
                Some(variants)
            }
            None => None,
        }
    } else {
        None
    };

    let fa_vcf_path = args.output.join("full_alignment.vcf");
    let mut fa_buf: Vec<u8> = Vec::new();
    let fa_inputs = FaPipelineInputs {
        bam_path: &bam,
        fasta_path: &fasta,
        region: Region {
            contig: contig.clone(),
            start: ctg_start.saturating_sub(1),
            end: ctg_end,
        },
        candidates: candidates_0b,
        phased_variants,
        server_endpoint: args.server.clone(),
        model_id: args.fa_model.clone(),
        sample_name: args.sample_name.clone(),
        clair3_version: args.clair3_version.clone(),
        matrix_depth,
        max_indel_length: 50,
        min_mq: 5,
        output_config: OutputWithConfig {
            is_show_reference: false,
            pileup: false,
            add_indel_length: true,
            ..OutputWithConfig::default()
        },
        batch_size: 32,
        aead_key: aead_key.clone(),
        tls: tls.clone(),
        low_mem: args.low_mem,
        fa_window_bp: args.fa_window_bp,
        streaming: false,
    };
    run_fa_pipeline(&fa_inputs, &mut fa_buf)
        .await
        .context("fa pipeline")?;
    let fa_text = String::from_utf8(fa_buf).context("fa utf8")?;
    let mut f = BufWriter::new(File::create(&fa_vcf_path)?);
    f.write_all(fa_text.as_bytes())?;
    f.flush()?;
    tracing::info!("[3/4] wrote {}", fa_vcf_path.display());

    tracing::info!("[4/4] merge_vcf");
    let merge_vcf_path = args.output.join("merge_output.vcf");
    let mut merge_buf: Vec<u8> = Vec::new();
    let stats = merge_vcf(
        pileup_text.as_bytes(),
        fa_text.as_bytes(),
        &mut merge_buf,
        Some(&contig),
        MergeConfig {
            quality_score_for_pass: Some(2.0),
            print_ref: false,
        },
    )
    .context("merge_vcf")?;

    let merge_buf = if region_filter.is_unrestricted() {
        merge_buf
    } else {
        filter_vcf_by_region(&merge_buf, &region_filter)
    };

    let mut f = BufWriter::new(File::create(&merge_vcf_path)?);
    f.write_all(&merge_buf)?;
    f.flush()?;
    tracing::info!(
        "[4/4] wrote {} ({} FA + {} pileup = {} rows pre-region-filter; {} pileup overridden)",
        merge_vcf_path.display(),
        stats.fa_rows,
        stats.pileup_rows,
        stats.fa_rows + stats.pileup_rows,
        stats.pileup_overridden,
    );
    println!(
        "Wrote {} ({} variants).",
        merge_vcf_path.display(),
        stats.fa_rows + stats.pileup_rows,
    );
    Ok(())
}
