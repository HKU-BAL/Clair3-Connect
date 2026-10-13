use anyhow::{anyhow, Context, Result};
use rust_htslib::bam::{self, IndexedReader, Read as _};
use rust_htslib::faidx;
use serde::Serialize;
use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::path::Path;

pub const EVIDENCE_FILTER_FLAGS: u16 = 4 | 8 | 256 | 512 | 1024 | 2048;

const BLOCK_SPAN: i64 = 50_000;

const MAX_BASE_QUAL: u8 = 40;

#[derive(Debug, Clone, Serialize)]
pub struct SiteQuery {
    pub contig: String,
    pub pos: i64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ref_allele: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub alt_allele: Option<String>,
}

impl SiteQuery {
    pub fn parse(s: &str) -> Result<Self> {
        let parts: Vec<&str> = s.trim().split(':').collect();
        let (contig, pos, r, a) = match parts.as_slice() {
            [c, p] => (*c, *p, None, None),
            [c, p, r, a] => (*c, *p, Some(r.to_ascii_uppercase()), Some(a.to_ascii_uppercase())),
            _ => return Err(anyhow!("site {s:?}: expected chrom:pos or chrom:pos:REF:ALT")),
        };
        let pos: i64 = pos
            .replace(',', "")
            .parse()
            .with_context(|| format!("site {s:?}: bad position"))?;
        if pos < 1 {
            return Err(anyhow!("site {s:?}: position must be >= 1"));
        }
        Ok(Self { contig: contig.to_string(), pos, ref_allele: r, alt_allele: a })
    }

    pub fn label(&self) -> String {
        format!("{}:{}", self.contig, self.pos)
    }
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct EvidenceParams {
    pub min_mapq: u8,
    pub min_baseq: u8,
    pub indel_qual: u8,
    pub threads: usize,
}

impl Default for EvidenceParams {
    fn default() -> Self {
        Self { min_mapq: 5, min_baseq: 0, indel_qual: 20, threads: 4 }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Obs {
    Base(u8),
    Ins,
    Del,
    Spanned,
}

impl Obs {
    fn label(&self) -> String {
        match self {
            Obs::Base(b) => (*b as char).to_string(),
            Obs::Ins => "INS".into(),
            Obs::Del => "DEL".into(),
            Obs::Spanned => "*".into(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct FragObs {
    obs: Obs,
    qual: u8,
    conflict: bool,
}

#[derive(Debug, Default, Clone)]
pub struct SiteObservations {
    pub raw_reads: u32,
    pub low_mapq_reads: u32,
    pub low_baseq_obs: u32,
    mapq_sum: u64,
    frags: HashMap<u64, FragObs>,
}

impl SiteObservations {
    fn add(&mut self, frag: u64, obs: Obs, qual: u8) {
        self.frags
            .entry(frag)
            .and_modify(|f| {
                if f.obs == obs {
                    f.qual = f.qual.max(qual);
                } else {
                    f.conflict = true;
                }
            })
            .or_insert(FragObs { obs, qual, conflict: false });
    }

    pub fn usable_depth(&self) -> usize {
        self.frags.values().filter(|f| !f.conflict).count()
    }

    fn usable(&self) -> impl Iterator<Item = (&u64, &FragObs)> {
        self.frags.iter().filter(|(_, f)| !f.conflict)
    }
}

fn frag_key(qname: &[u8]) -> u64 {
    let mut h = DefaultHasher::new();
    qname.hash(&mut h);
    h.finish()
}

fn err_rate(q: u8) -> f64 {
    10f64.powf(-(q.min(MAX_BASE_QUAL).max(1) as f64) / 10.0)
}

pub fn collect_observations(
    bam_path: &Path,
    sites: &[SiteQuery],
    params: &EvidenceParams,
) -> Result<Vec<SiteObservations>> {
    let mut order: Vec<usize> = (0..sites.len()).collect();
    order.sort_by(|&a, &b| (&sites[a].contig, sites[a].pos).cmp(&(&sites[b].contig, sites[b].pos)));
    let mut blocks: Vec<Vec<usize>> = Vec::new();
    for idx in order {
        match blocks.last_mut() {
            Some(b)
                if sites[b[0]].contig == sites[idx].contig
                    && sites[idx].pos - sites[b[0]].pos <= BLOCK_SPAN =>
            {
                b.push(idx)
            }
            _ => blocks.push(vec![idx]),
        }
    }

    let n_threads = params.threads.max(1).min(blocks.len().max(1));
    let chunks: Vec<Vec<Vec<usize>>> = {
        let mut c: Vec<Vec<Vec<usize>>> = vec![Vec::new(); n_threads];
        for (i, b) in blocks.into_iter().enumerate() {
            c[i % n_threads].push(b);
        }
        c
    };

    let results: Vec<Result<Vec<(usize, SiteObservations)>>> = std::thread::scope(|s| {
        let handles: Vec<_> = chunks
            .into_iter()
            .map(|chunk| {
                s.spawn(move || -> Result<Vec<(usize, SiteObservations)>> {
                    let mut reader = IndexedReader::from_path(bam_path)
                        .with_context(|| format!("open BAM {}", bam_path.display()))?;
                    let mut out = Vec::new();
                    for block in chunk {
                        out.extend(process_block(&mut reader, sites, &block, params)?);
                    }
                    Ok(out)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().expect("evidence worker panicked")).collect()
    });

    let mut per_site: Vec<Option<SiteObservations>> = vec![None; sites.len()];
    for r in results {
        for (i, o) in r? {
            per_site[i] = Some(o);
        }
    }
    Ok(per_site.into_iter().map(|o| o.unwrap_or_default()).collect())
}

fn process_block(
    reader: &mut IndexedReader,
    sites: &[SiteQuery],
    block: &[usize],
    params: &EvidenceParams,
) -> Result<Vec<(usize, SiteObservations)>> {
    let contig = &sites[block[0]].contig;
    let pos0: Vec<i64> = block.iter().map(|&i| sites[i].pos - 1).collect();
    let mut obs: Vec<SiteObservations> = vec![SiteObservations::default(); block.len()];

    let tid = match reader.header().tid(contig.as_bytes()) {
        Some(t) => t,
        None => return Ok(block.iter().copied().zip(obs).collect()),
    };
    let (start, end) = (pos0[0], *pos0.last().unwrap() + 1);
    reader
        .fetch((tid, start, end))
        .with_context(|| format!("fetch {contig}:{start}-{end}"))?;

    let mut record = bam::Record::new();
    while let Some(r) = reader.read(&mut record) {
        r.context("read BAM record")?;
        if record.flags() & EVIDENCE_FILTER_FLAGS != 0 {
            continue;
        }
        let rstart = record.pos();
        let cigar = record.cigar();
        let rend = cigar.end_pos();
        let first = pos0.partition_point(|&p| p < rstart);
        if first >= pos0.len() || pos0[first] >= rend {
            continue;
        }
        let mapq = record.mapq();
        let low_mq = mapq < params.min_mapq;
        let frag = frag_key(record.qname());
        let seq = record.seq();
        let qual = record.qual();

        let mut k = first;
        let mut rpos = rstart;
        let mut qpos: usize = 0;
        let ops: Vec<_> = cigar.iter().collect();
        for (oi, op) in ops.iter().enumerate() {
            if k >= pos0.len() {
                break;
            }
            let len = op.len() as i64;
            match op.char() {
                'M' | '=' | 'X' => {
                    let block_end = rpos + len;
                    while k < pos0.len() && pos0[k] < block_end {
                        let p = pos0[k];
                        if p >= rpos {
                            let o = &mut obs[k];
                            o.raw_reads += 1;
                            o.mapq_sum += mapq as u64;
                            if low_mq {
                                o.low_mapq_reads += 1;
                            } else {
                                let qi = qpos + (p - rpos) as usize;
                                let bq = qual[qi];
                                let next = if p == block_end - 1 { ops.get(oi + 1).map(|c| c.char()) } else { None };
                                match next {
                                    Some('I') => o.add(frag, Obs::Ins, params.indel_qual),
                                    Some('D') => o.add(frag, Obs::Del, params.indel_qual),
                                    _ => {
                                        if bq < params.min_baseq {
                                            o.low_baseq_obs += 1;
                                        } else {
                                            let b = seq[qi].to_ascii_uppercase();
                                            o.add(frag, Obs::Base(b), bq);
                                        }
                                    }
                                }
                            }
                        }
                        k += 1;
                    }
                    rpos = block_end;
                    qpos += len as usize;
                }
                'D' => {
                    let block_end = rpos + len;
                    while k < pos0.len() && pos0[k] < block_end {
                        if pos0[k] >= rpos {
                            let o = &mut obs[k];
                            o.raw_reads += 1;
                            o.mapq_sum += mapq as u64;
                            if low_mq {
                                o.low_mapq_reads += 1;
                            } else {
                                o.add(frag, Obs::Spanned, params.indel_qual);
                            }
                        }
                        k += 1;
                    }
                    rpos = block_end;
                }
                'N' => {
                    let block_end = rpos + len;
                    while k < pos0.len() && pos0[k] < block_end {
                        k += 1;
                    }
                    rpos = block_end;
                }
                'I' | 'S' => qpos += len as usize,
                _ => {}
            }
        }
    }
    Ok(block.iter().copied().zip(obs).collect())
}

#[derive(Debug, Default, Clone, Serialize)]
pub struct AlleleCounts {
    #[allow(missing_docs)]
    pub a: u32,
    #[allow(missing_docs)]
    pub c: u32,
    #[allow(missing_docs)]
    pub g: u32,
    #[allow(missing_docs)]
    pub t: u32,
    #[allow(missing_docs)]
    pub n: u32,
    pub ins: u32,
    pub del: u32,
    pub spanned: u32,
}

#[derive(Debug, Clone, Serialize)]
pub struct SiteEvidence {
    pub site: String,
    pub ref_base: String,
    pub raw_depth: u32,
    pub low_mapq_reads: u32,
    pub mean_mapq: f64,
    pub depth: u32,
    pub counts: AlleleCounts,
    pub top_non_ref: String,
    pub top_non_ref_fraction: f64,
}

pub fn summarize(site: &SiteQuery, ref_base: u8, o: &SiteObservations) -> SiteEvidence {
    let ref_base = ref_base.to_ascii_uppercase();
    let mut counts = AlleleCounts::default();
    let mut tally: HashMap<Obs, u32> = HashMap::new();
    let mut depth = 0u32;
    for (_, f) in o.usable() {
        depth += 1;
        match f.obs {
            Obs::Base(b'A') => counts.a += 1,
            Obs::Base(b'C') => counts.c += 1,
            Obs::Base(b'G') => counts.g += 1,
            Obs::Base(b'T') => counts.t += 1,
            Obs::Base(_) => counts.n += 1,
            Obs::Ins => counts.ins += 1,
            Obs::Del => counts.del += 1,
            Obs::Spanned => counts.spanned += 1,
        }
        if f.obs != Obs::Base(ref_base) && f.obs != Obs::Base(b'N') {
            *tally.entry(f.obs).or_default() += 1;
        }
    }
    let top = tally.iter().max_by_key(|(o, c)| (**c, o.label())).map(|(o, c)| (o.label(), *c));
    let (top_non_ref, top_count) = top.unwrap_or((".".into(), 0));
    SiteEvidence {
        site: site.label(),
        ref_base: (ref_base as char).to_string(),
        raw_depth: o.raw_reads,
        low_mapq_reads: o.low_mapq_reads,
        mean_mapq: if o.raw_reads > 0 {
            ((o.mapq_sum as f64 / o.raw_reads as f64) * 10.0).round() / 10.0
        } else {
            0.0
        },
        depth,
        counts,
        top_non_ref,
        top_non_ref_fraction: if depth > 0 {
            ((top_count as f64 / depth as f64) * 10_000.0).round() / 10_000.0
        } else {
            0.0
        },
    }
}

pub fn reference_bases(fasta: &faidx::Reader, sites: &[SiteQuery]) -> Vec<u8> {
    sites
        .iter()
        .map(|s| {
            let p = (s.pos - 1) as usize;
            fasta
                .fetch_seq(&s.contig, p, p)
                .ok()
                .and_then(|v| v.first().copied())
                .map(|b| b.to_ascii_uppercase())
                .unwrap_or(b'N')
        })
        .collect()
}

#[derive(Debug, Clone, Serialize)]
pub struct PairEvidence {
    pub a: String,
    pub b: String,
    pub cis: u32,
    pub trans: u32,
    pub other: u32,
    pub llr10: f64,
    pub pq: u32,
    pub relation: &'static str,
}

#[derive(Debug, Clone, Serialize)]
pub struct PhaseBlock {
    pub sites: Vec<String>,
    pub hap1: Vec<String>,
    pub hap2: Vec<String>,
    pub min_pq: u32,
}

#[derive(Debug, Clone, Serialize)]
pub struct PhaseSiteSupport {
    pub site: String,
    #[allow(missing_docs)]
    #[serde(rename = "ref")]
    pub ref_allele: String,
    #[allow(missing_docs)]
    pub alt: String,
    pub ref_count: u32,
    pub alt_count: u32,
    pub other_count: u32,
    pub heterozygous: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub warning: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PhaseResult {
    pub status: &'static str,
    pub sites: Vec<PhaseSiteSupport>,
    pub pairs: Vec<PairEvidence>,
    pub blocks: Vec<PhaseBlock>,
    pub excluded_sites: Vec<String>,
    pub min_pq: u32,
}

pub fn phase_sites(
    bam_path: &Path,
    sites: &[SiteQuery],
    ref_bases: &[u8],
    alt_bases: &[u8],
    params: &EvidenceParams,
    min_pq: u32,
) -> Result<PhaseResult> {
    let obs = collect_observations(bam_path, sites, params)?;
    let n = sites.len();

    let mut frags: HashMap<u64, Vec<Option<(u8, f64)>>> = HashMap::new();
    let mut support = Vec::with_capacity(n);
    for (i, o) in obs.iter().enumerate() {
        let (mut r, mut a, mut x) = (0u32, 0u32, 0u32);
        for (k, f) in o.usable() {
            let code = match f.obs {
                Obs::Base(b) if b == ref_bases[i] => 0u8,
                Obs::Base(b) if b == alt_bases[i] => 1u8,
                _ => 2u8,
            };
            match code {
                0 => r += 1,
                1 => a += 1,
                _ => x += 1,
            }
            frags.entry(*k).or_insert_with(|| vec![None; n])[i] = Some((code, err_rate(f.qual)));
        }
        let tot = r + a;
        let af = if tot > 0 { a as f64 / tot as f64 } else { 0.0 };
        let heterozygous = tot > 0 && (0.2..=0.8).contains(&af);
        let warning = if tot == 0 {
            Some("no fragment shows REF or ALT".to_string())
        } else {
            let af = a as f64 / tot as f64;
            if af < 0.2 {
                Some(format!("ALT fraction {af:.2}: site does not look heterozygous (mostly REF)"))
            } else if af > 0.8 {
                Some(format!("ALT fraction {af:.2}: site does not look heterozygous (mostly ALT)"))
            } else {
                None
            }
        };
        support.push(PhaseSiteSupport {
            site: sites[i].label(),
            ref_allele: (ref_bases[i] as char).to_string(),
            alt: (alt_bases[i] as char).to_string(),
            ref_count: r,
            alt_count: a,
            other_count: x,
            heterozygous,
            warning,
        });
    }

    let het: Vec<bool> = support.iter().map(|s| s.heterozygous).collect();
    let mut llr = vec![vec![0f64; n]; n];
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if !(het[i] && het[j]) {
                continue;
            }
            let (mut cis, mut trans, mut other, mut l) = (0u32, 0u32, 0u32, 0f64);
            for v in frags.values() {
                if let (Some((ci, ei)), Some((cj, ej))) = (v[i], v[j]) {
                    if ci == 2 || cj == 2 {
                        other += 1;
                        continue;
                    }
                    let q = (ei + ej - 2.0 * ei * ej).clamp(0.01, 0.49);
                    let w = ((1.0 - q) / q).log10();
                    if ci == cj {
                        cis += 1;
                        l += w;
                    } else {
                        trans += 1;
                        l -= w;
                    }
                }
            }
            if cis + trans + other == 0 {
                continue;
            }
            let pq = pq_from_llr(l);
            let relation = if pq < min_pq { "undetermined" } else if l > 0.0 { "cis" } else { "trans" };
            llr[i][j] = l;
            llr[j][i] = l;
            pairs.push(PairEvidence {
                a: sites[i].label(),
                b: sites[j].label(),
                cis,
                trans,
                other,
                llr10: (l * 100.0).round() / 100.0,
                pq,
                relation,
            });
        }
    }

    let mut assigned: Vec<Option<u8>> = vec![None; n];
    let mut blocks = Vec::new();
    for root in 0..n {
        if assigned[root].is_some() || !het[root] {
            continue;
        }
        assigned[root] = Some(1);
        let mut members = vec![root];
        let mut weakest = 99u32;
        loop {
            let mut best: Option<(usize, usize, f64)> = None;
            for &u in &members {
                for v in 0..n {
                    if assigned[v].is_some() || pq_from_llr(llr[u][v]) < min_pq {
                        continue;
                    }
                    if best.map(|b| llr[u][v].abs() > b.2.abs()).unwrap_or(true) {
                        best = Some((u, v, llr[u][v]));
                    }
                }
            }
            match best {
                Some((u, v, l)) => {
                    let hu = assigned[u].unwrap();
                    assigned[v] = Some(if l > 0.0 { hu } else { 1 - hu });
                    weakest = weakest.min(pq_from_llr(l));
                    members.push(v);
                }
                None => break,
            }
        }
        members.sort();
        let allele = |i: usize, h: u8| -> String {
            let b = if h == 1 { alt_bases[i] } else { ref_bases[i] };
            (b as char).to_string()
        };
        blocks.push(PhaseBlock {
            sites: members.iter().map(|&i| sites[i].label()).collect(),
            hap1: members.iter().map(|&i| allele(i, assigned[i].unwrap())).collect(),
            hap2: members.iter().map(|&i| allele(i, 1 - assigned[i].unwrap())).collect(),
            min_pq: weakest,
        });
    }

    let excluded_sites: Vec<String> = (0..n).filter(|&i| !het[i]).map(|i| sites[i].label()).collect();
    let n_het = het.iter().filter(|&&h| h).count();
    let status = if n_het >= 2 && blocks.len() == 1 {
        "phased"
    } else if blocks.iter().any(|b| b.sites.len() > 1) {
        "partially_phased"
    } else {
        "unphased"
    };
    Ok(PhaseResult { status, sites: support, pairs, blocks, excluded_sites, min_pq })
}

fn pq_from_llr(l: f64) -> u32 {
    let a = l.abs();
    let pq = 10.0 * (1.0 + 10f64.powf(a)).log10();
    (pq.round() as u32).min(99)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_sites() {
        let s = SiteQuery::parse("chr19:44,908,684:T:C").unwrap();
        assert_eq!((s.contig.as_str(), s.pos), ("chr19", 44908684));
        assert_eq!(s.alt_allele.as_deref(), Some("C"));
        assert!(SiteQuery::parse("chr19").is_err());
        assert!(SiteQuery::parse("chr19:0").is_err());
    }

    #[test]
    fn counts_and_mate_merging() {
        let site = SiteQuery::parse("chr1:100").unwrap();
        let mut o = SiteObservations { raw_reads: 5, ..Default::default() };
        o.add(1, Obs::Base(b'A'), 30);
        o.add(1, Obs::Base(b'A'), 35);
        o.add(2, Obs::Base(b'G'), 30);
        o.add(3, Obs::Base(b'G'), 30);
        o.add(3, Obs::Base(b'A'), 30);
        let e = summarize(&site, b'A', &o);
        assert_eq!((e.depth, e.counts.a, e.counts.g), (2, 1, 1));
        assert_eq!((e.top_non_ref.as_str(), e.top_non_ref_fraction), ("G", 0.5));
    }

    #[test]
    fn pq_is_monotone() {
        assert_eq!(pq_from_llr(0.0), 3);
        assert!(pq_from_llr(2.0) >= 20);
        assert_eq!(pq_from_llr(50.0), 99);
    }
}
