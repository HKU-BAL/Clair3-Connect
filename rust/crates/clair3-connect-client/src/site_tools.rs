use anyhow::{anyhow, bail, Context, Result};
use clair3_connect_core::callability::{scan_callability, CallabilitySummary, CallableParams};
use clair3_connect_core::site_evidence::{
    collect_observations, phase_sites, reference_bases, summarize, EvidenceParams, SiteEvidence, SiteQuery,
};
use rust_htslib::faidx;
use serde::Serialize;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

const SITES_USAGE: &str = "clair3-connect sites --bam BAM --ref REF.fa \\
        (--sites chr:pos[:REF:ALT][,...] | --sites-file FILE) [--vcf CALLS.vcf[.gz]] \\
        [--min-mapq 5] [--min-baseq 0] [--threads 4] [--format json|tsv] [--output FILE]

Descriptive read evidence at each position: reads covering it, reads below
--min-mapq, mean MAPQ, usable fragments (overlapping mates merged), and the
base distribution (A/C/G/T/N, insertion or deletion anchored at the position,
deletion spanning it). With --vcf, Clair3's record at the position (if any) is
reported alongside. No genotype is inferred.";

const CALLABLE_USAGE: &str = "clair3-connect callable --bam BAM \\
        (--region chr:start-end | --regions chr:s-e[,...] | --bed FILE) --output-dir DIR \\
        [--platform ont|hifi|ilmn] [--min-depth D] [--min-mapq 5] [--poor-mapq 20] [--max-poor-mapq-fraction 0.1] \\
        [--max-depth-multiple 2.5] [--threads 4]

Writes DIR/callable.bed (merged intervals: chrom, start, end, state) and
DIR/callability.json (bases per state, callable fraction, median depth).
States: NO_COVERAGE (no read), POOR_MAPPING_QUALITY (>= fraction of covering
reads below --poor-mapq), LOW_COVERAGE (< --min-depth reads with MAPQ >=
--min-mapq; default 8 for ont/hifi, 10 for ilmn), EXCESSIVE_COVERAGE (>
multiple x the 1-Mb window median), CALLABLE.";

const PHASE_USAGE: &str = "clair3-connect phase --bam BAM --ref REF.fa \\
        --sites chr:pos[:REF:ALT][,...] [--vcf CALLS.vcf[.gz]] \\
        [--min-mapq 5] [--min-baseq 0] [--min-pq 20] [--threads 4] [--output FILE]

Phases SNV sites from fragments spanning two or more of them (mates merged).
ALT alleles come from the site spec, else from Clair3's record at the site
(--vcf), else from the most frequent non-reference base in the reads.";

#[derive(Default)]
struct Opts {
    bam: Option<PathBuf>,
    fasta: Option<PathBuf>,
    sites: Vec<SiteQuery>,
    vcf: Option<PathBuf>,
    regions: Vec<(String, u64, u64)>,
    output: Option<PathBuf>,
    output_dir: Option<PathBuf>,
    format: String,
    ev: EvidenceParams,
    cp: CallableParams,
    min_pq: u32,
    min_depth: Option<u32>,
    platform: Option<clair3_connect_core::platform::Platform>,
}

fn parse(argv: &[String], usage: &str) -> Result<Opts> {
    let mut o = Opts { format: "json".into(), min_pq: 20, ..Default::default() };
    let mut it = argv.iter();
    while let Some(a) = it.next() {
        let mut val = || it.next().cloned().ok_or_else(|| anyhow!("flag {a} needs a value"));
        match a.as_str() {
            "--bam" => o.bam = Some(PathBuf::from(val()?)),
            "--ref" => o.fasta = Some(PathBuf::from(val()?)),
            "--sites" => {
                for s in val()?.split(',').filter(|s| !s.trim().is_empty()) {
                    o.sites.push(SiteQuery::parse(s)?);
                }
            }
            "--sites-file" => {
                let f = std::fs::File::open(val()?).context("open --sites-file")?;
                for line in BufReader::new(f).lines() {
                    let line = line?;
                    let t = line.split('\t').next().unwrap_or("").trim().to_string();
                    if !t.is_empty() && !t.starts_with('#') {
                        o.sites.push(SiteQuery::parse(&t)?);
                    }
                }
            }
            "--vcf" => o.vcf = Some(PathBuf::from(val()?)),
            "--region" | "--regions" => {
                for r in val()?.split(',').filter(|s| !s.trim().is_empty()) {
                    let (c, s, e) = parse_region(r)?;
                    o.regions.push((c, (s - 1) as u64, e as u64));
                }
            }
            "--bed" => {
                let f = std::fs::File::open(val()?).context("open --bed")?;
                for line in BufReader::new(f).lines() {
                    let line = line?;
                    let c: Vec<&str> = line.split('\t').collect();
                    if c.len() >= 3 && !line.starts_with('#') {
                        o.regions.push((c[0].to_string(), c[1].parse()?, c[2].parse()?));
                    }
                }
            }
            "--output" => o.output = Some(PathBuf::from(val()?)),
            "--output-dir" => o.output_dir = Some(PathBuf::from(val()?)),
            "--format" => o.format = val()?,
            "--min-mapq" => {
                let m: u8 = val()?.parse()?;
                o.ev.min_mapq = m;
                o.cp.min_mapq = m;
            }
            "--min-baseq" => o.ev.min_baseq = val()?.parse()?,
            "--threads" => {
                let t: usize = val()?.parse()?;
                o.ev.threads = t;
                o.cp.threads = t;
            }
            "--min-depth" => o.min_depth = Some(val()?.parse()?),
            "--platform" => o.platform = Some(clair3_connect_core::platform::Platform::from_str_or_ont(&val()?)),
            "--poor-mapq" => o.cp.poor_mapq = val()?.parse()?,
            "--max-poor-mapq-fraction" => o.cp.max_poor_mapq_fraction = val()?.parse()?,
            "--max-depth-multiple" => o.cp.max_depth_multiple = val()?.parse()?,
            "--min-pq" => o.min_pq = val()?.parse()?,
            "-h" | "--help" => {
                eprintln!("{usage}");
                std::process::exit(0);
            }
            other => bail!("unknown flag {other}\n\n{usage}"),
        }
    }
    if let Some(p) = o.platform {
        let d = CallableParams::for_platform(p);
        o.cp.min_depth = d.min_depth;
    }
    if let Some(d) = o.min_depth {
        o.cp.min_depth = d;
    }
    Ok(o)
}

pub fn parse_region(s: &str) -> Result<(String, i64, i64)> {
    let s = s.trim().replace(',', "");
    let (c, se) = s.rsplit_once(':').ok_or_else(|| anyhow!("region {s:?}: expected chr:start-end"))?;
    let (a, b) = se.split_once('-').ok_or_else(|| anyhow!("region {s:?}: expected chr:start-end"))?;
    Ok((c.to_string(), a.parse()?, b.parse()?))
}

fn need<T>(v: Option<T>, what: &str, usage: &str) -> Result<T> {
    v.ok_or_else(|| anyhow!("{what} is required\n\n{usage}"))
}

fn open_text(path: &Path) -> Result<Box<dyn BufRead>> {
    let f = std::fs::File::open(path).with_context(|| format!("open {}", path.display()))?;
    if path.extension().map(|e| e == "gz").unwrap_or(false) {
        Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(f))))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

#[derive(Debug, Clone, Serialize)]
struct VcfRecord {
    pos: i64,
    #[serde(rename = "ref")]
    ref_allele: String,
    alt: String,
    filter: String,
    gt: String,
    qual: f64,
}

fn vcf_records_at(path: &Path, sites: &[SiteQuery]) -> Result<Vec<Option<VcfRecord>>> {
    let mut want: HashMap<(String, i64), Vec<usize>> = HashMap::new();
    for (i, s) in sites.iter().enumerate() {
        want.entry((s.contig.clone(), s.pos)).or_default().push(i);
    }
    let mut out = vec![None; sites.len()];
    for line in open_text(path)?.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 10 {
            continue;
        }
        let pos: i64 = c[1].parse().unwrap_or(0);
        let Some(idx) = want.get(&(c[0].to_string(), pos)) else { continue };
        let fmt: Vec<&str> = c[8].split(':').collect();
        let val: Vec<&str> = c[9].split(':').collect();
        let gt = fmt
            .iter()
            .position(|f| *f == "GT")
            .and_then(|i| val.get(i))
            .map(|g| g.replace('|', "/"))
            .unwrap_or_else(|| "./.".into());
        let rec = VcfRecord {
            pos,
            ref_allele: c[3].into(),
            alt: c[4].into(),
            filter: c[6].into(),
            gt,
            qual: c[5].parse().unwrap_or(0.0),
        };
        for &i in idx {
            out[i] = Some(rec.clone());
        }
    }
    Ok(out)
}

fn write_out(o: &Opts, text: &str) -> Result<()> {
    match &o.output {
        Some(p) => std::fs::write(p, text).with_context(|| format!("write {}", p.display())),
        None => {
            std::io::stdout().lock().write_all(text.as_bytes())?;
            Ok(())
        }
    }
}

#[derive(Serialize)]
struct SiteReport {
    #[serde(flatten)]
    query: SiteQuery,
    clair3_record: Option<VcfRecord>,
    evidence: SiteEvidence,
}

pub fn run_sites(argv: &[String]) -> Result<()> {
    let o = parse(argv, SITES_USAGE)?;
    let bam = need(o.bam.clone(), "--bam", SITES_USAGE)?;
    let fasta_path = need(o.fasta.clone(), "--ref", SITES_USAGE)?;
    if o.sites.is_empty() {
        bail!("no sites given (--sites or --sites-file)\n\n{SITES_USAGE}");
    }
    let fasta = faidx::Reader::from_path(&fasta_path).with_context(|| format!("open {}", fasta_path.display()))?;
    let refs = reference_bases(&fasta, &o.sites);
    let obs = collect_observations(&bam, &o.sites, &o.ev)?;
    let recs = match &o.vcf {
        Some(v) => vcf_records_at(v, &o.sites)?,
        None => vec![None; o.sites.len()],
    };
    let reports: Vec<SiteReport> = o
        .sites
        .iter()
        .enumerate()
        .map(|(i, s)| SiteReport { query: s.clone(), clair3_record: recs[i].clone(), evidence: summarize(s, refs[i], &obs[i]) })
        .collect();

    if o.format == "tsv" {
        let mut t = String::from(
            "site\tclair3_filter\tclair3_gt\tref_base\traw_depth\tlow_mapq_reads\tmean_mapq\tdepth\tA\tC\tG\tT\tN\tINS\tDEL\tSPANNED\ttop_non_ref\ttop_non_ref_fraction\n",
        );
        for r in &reports {
            let e = &r.evidence;
            let k = &e.counts;
            t.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                e.site,
                r.clair3_record.as_ref().map(|v| v.filter.as_str()).unwrap_or("."),
                r.clair3_record.as_ref().map(|v| v.gt.as_str()).unwrap_or("."),
                e.ref_base,
                e.raw_depth,
                e.low_mapq_reads,
                e.mean_mapq,
                e.depth,
                k.a,
                k.c,
                k.g,
                k.t,
                k.n,
                k.ins,
                k.del,
                k.spanned,
                e.top_non_ref,
                e.top_non_ref_fraction
            ));
        }
        return write_out(&o, &t);
    }
    #[derive(Serialize)]
    struct Out<'a> {
        tool: &'static str,
        version: &'static str,
        evidence_params: EvidenceParams,
        note: &'static str,
        sites: &'a [SiteReport],
    }
    let out = Out {
        tool: "clair3-connect sites",
        version: env!("CARGO_PKG_VERSION"),
        evidence_params: o.ev,
        note: "descriptive read evidence; genotypes come from Clair3 (clair3_record)",
        sites: &reports,
    };
    write_out(&o, &(serde_json::to_string_pretty(&out)? + "\n"))
}

pub fn write_callability(
    bam: &Path,
    regions: &[(String, u64, u64)],
    params: &CallableParams,
    out_dir: &Path,
) -> Result<CallabilitySummary> {
    let t0 = std::time::Instant::now();
    let (ivs, summary) = scan_callability(bam, regions, params)?;
    std::fs::create_dir_all(out_dir)?;
    let mut f = BufWriter::new(std::fs::File::create(out_dir.join("callable.bed"))?);
    for iv in &ivs {
        writeln!(f, "{}\t{}\t{}\t{}", iv.contig, iv.start, iv.end, iv.state.label())?;
    }
    f.flush()?;
    let mut j = serde_json::to_value(&summary)?;
    j["regions"] = serde_json::json!(regions.iter().map(|(c, s, e)| format!("{c}:{}-{e}", s + 1)).collect::<Vec<_>>());
    j["elapsed_s"] = serde_json::json!((t0.elapsed().as_secs_f64() * 100.0).round() / 100.0);
    std::fs::write(out_dir.join("callability.json"), serde_json::to_string_pretty(&j)? + "\n")?;
    Ok(summary)
}

pub fn run_callable(argv: &[String]) -> Result<()> {
    let o = parse(argv, CALLABLE_USAGE)?;
    let bam = need(o.bam.clone(), "--bam", CALLABLE_USAGE)?;
    let dir = need(o.output_dir.clone(), "--output-dir", CALLABLE_USAGE)?;
    if o.regions.is_empty() {
        bail!("no regions given\n\n{CALLABLE_USAGE}");
    }
    let s = write_callability(&bam, &o.regions, &o.cp, &dir)?;
    println!("{}", serde_json::to_string(&s)?);
    Ok(())
}

pub fn run_phase(argv: &[String]) -> Result<()> {
    let o = parse(argv, PHASE_USAGE)?;
    let bam = need(o.bam.clone(), "--bam", PHASE_USAGE)?;
    let fasta_path = need(o.fasta.clone(), "--ref", PHASE_USAGE)?;
    if o.sites.len() < 2 {
        bail!("phase needs at least two sites\n\n{PHASE_USAGE}");
    }
    let fasta = faidx::Reader::from_path(&fasta_path).with_context(|| format!("open {}", fasta_path.display()))?;
    let refs = reference_bases(&fasta, &o.sites);
    let recs = match &o.vcf {
        Some(v) => vcf_records_at(v, &o.sites)?,
        None => vec![None; o.sites.len()],
    };
    let mut alts = Vec::with_capacity(o.sites.len());
    let mut need_reads = false;
    for (i, s) in o.sites.iter().enumerate() {
        if let Some(r) = &s.ref_allele {
            if r.len() != 1 || r.as_bytes()[0] != refs[i] {
                bail!("site {}: REF {r} does not match the reference base {}", s.label(), refs[i] as char);
            }
        }
        let alt = s
            .alt_allele
            .clone()
            .or_else(|| recs[i].as_ref().map(|r| r.alt.split(',').next().unwrap_or("").to_string()));
        match alt {
            Some(a) if a.len() == 1 => alts.push(a.as_bytes()[0]),
            Some(a) => bail!("site {}: only SNVs can be phased (ALT {a})", s.label()),
            None => {
                need_reads = true;
                alts.push(0);
            }
        }
    }
    if need_reads {
        let obs = collect_observations(&bam, &o.sites, &o.ev)?;
        for i in 0..o.sites.len() {
            if alts[i] == 0 {
                let e = summarize(&o.sites[i], refs[i], &obs[i]);
                if e.top_non_ref.len() != 1 {
                    bail!("site {}: no SNV alternate allele found in the reads", o.sites[i].label());
                }
                alts[i] = e.top_non_ref.as_bytes()[0];
            }
        }
    }
    let res = phase_sites(&bam, &o.sites, &refs, &alts, &o.ev, o.min_pq)?;
    write_out(&o, &(serde_json::to_string_pretty(&res)? + "\n"))
}
