use anyhow::{Context, Result};
use rust_htslib::bam::{self, IndexedReader, Read as _};
use serde::Serialize;
use std::path::Path;

pub const CALLABILITY_FILTER_FLAGS: u16 = 4 | 8 | 256 | 512 | 1024 | 2048;

const WINDOW: u64 = 1_000_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[allow(missing_docs)]
pub enum CallState {
    Callable,
    NoCoverage,
    LowCoverage,
    PoorMappingQuality,
    ExcessiveCoverage,
}

impl CallState {
    pub fn label(&self) -> &'static str {
        match self {
            CallState::Callable => "CALLABLE",
            CallState::NoCoverage => "NO_COVERAGE",
            CallState::LowCoverage => "LOW_COVERAGE",
            CallState::PoorMappingQuality => "POOR_MAPPING_QUALITY",
            CallState::ExcessiveCoverage => "EXCESSIVE_COVERAGE",
        }
    }
    const ALL: [CallState; 5] = [
        CallState::Callable,
        CallState::NoCoverage,
        CallState::LowCoverage,
        CallState::PoorMappingQuality,
        CallState::ExcessiveCoverage,
    ];
}

#[derive(Debug, Clone, Copy, Serialize)]
pub struct CallableParams {
    pub min_depth: u32,
    pub min_mapq: u8,
    pub poor_mapq: u8,
    pub max_poor_mapq_fraction: f64,
    pub max_depth_multiple: f64,
    pub threads: usize,
}

impl CallableParams {
    pub fn for_platform(platform: crate::platform::Platform) -> Self {
        let min_depth = match platform {
            crate::platform::Platform::Ilmn => 10,
            _ => 8,
        };
        Self { min_depth, ..Self::default() }
    }
}

impl Default for CallableParams {
    fn default() -> Self {
        Self {
            min_depth: 8,
            min_mapq: 5,
            poor_mapq: 20,
            max_poor_mapq_fraction: 0.1,
            max_depth_multiple: 2.5,
            threads: 4,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct CallableInterval {
    #[allow(missing_docs)]
    pub contig: String,
    #[allow(missing_docs)]
    pub start: u64,
    #[allow(missing_docs)]
    pub end: u64,
    #[allow(missing_docs)]
    pub state: CallState,
}

#[derive(Debug, Clone, Serialize)]
pub struct CallabilitySummary {
    pub params: CallableParams,
    pub total_bp: u64,
    pub bp_by_state: std::collections::BTreeMap<String, u64>,
    pub callable_fraction: f64,
    pub median_depth: f64,
}

pub fn scan_callability(
    bam_path: &Path,
    regions: &[(String, u64, u64)],
    p: &CallableParams,
) -> Result<(Vec<CallableInterval>, CallabilitySummary)> {
    let mut windows: Vec<(String, u64, u64)> = Vec::new();
    for (c, s, e) in regions {
        let mut a = *s;
        while a < *e {
            let b = (a + WINDOW).min(*e);
            windows.push((c.clone(), a, b));
            a = b;
        }
    }
    let n_threads = p.threads.max(1).min(windows.len().max(1));
    let results: Vec<Result<Vec<(usize, Vec<CallableInterval>, f64)>>> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..n_threads)
            .map(|t| {
                let windows = &windows;
                sc.spawn(move || -> Result<Vec<(usize, Vec<CallableInterval>, f64)>> {
                    let mut reader = IndexedReader::from_path(bam_path)
                        .with_context(|| format!("open BAM {}", bam_path.display()))?;
                    let mut out = Vec::new();
                    for i in (t..windows.len()).step_by(n_threads) {
                        let (iv, med) = scan_window(&mut reader, &windows[i], p)?;
                        out.push((i, iv, med));
                    }
                    Ok(out)
                })
            })
            .collect();
        hs.into_iter().map(|h| h.join().expect("callability worker panicked")).collect()
    });
    let mut per: Vec<(usize, Vec<CallableInterval>, f64)> = Vec::new();
    for r in results {
        per.extend(r?);
    }
    per.sort_by_key(|x| x.0);

    let mut merged: Vec<CallableInterval> = Vec::new();
    let mut medians: Vec<(f64, u64)> = Vec::new();
    for (i, ivs, med) in per {
        medians.push((med, windows[i].2 - windows[i].1));
        for iv in ivs {
            match merged.last_mut() {
                Some(m) if m.contig == iv.contig && m.end == iv.start && m.state == iv.state => m.end = iv.end,
                _ => merged.push(iv),
            }
        }
    }
    let mut bp_by_state = std::collections::BTreeMap::new();
    for s in CallState::ALL {
        bp_by_state.insert(s.label().to_string(), 0u64);
    }
    let mut total = 0u64;
    for iv in &merged {
        *bp_by_state.get_mut(iv.state.label()).unwrap() += iv.end - iv.start;
        total += iv.end - iv.start;
    }
    medians.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let half = medians.iter().map(|m| m.1).sum::<u64>() / 2;
    let mut acc = 0u64;
    let mut median_depth = 0.0;
    for (m, w) in &medians {
        acc += w;
        if acc > half {
            median_depth = *m;
            break;
        }
    }
    let callable = bp_by_state["CALLABLE"];
    Ok((
        merged,
        CallabilitySummary {
            params: *p,
            total_bp: total,
            bp_by_state,
            callable_fraction: if total > 0 { (callable as f64 / total as f64 * 1e4).round() / 1e4 } else { 0.0 },
            median_depth,
        },
    ))
}

fn scan_window(
    reader: &mut IndexedReader,
    w: &(String, u64, u64),
    p: &CallableParams,
) -> Result<(Vec<CallableInterval>, f64)> {
    let (contig, ws, we) = (w.0.as_str(), w.1 as i64, w.2 as i64);
    let len = (we - ws) as usize;
    let mut d_all = vec![0i32; len + 1];
    let mut d_use = vec![0i32; len + 1];
    let mut d_poor = vec![0i32; len + 1];
    if let Some(tid) = reader.header().tid(contig.as_bytes()) {
        reader.fetch((tid, ws, we)).with_context(|| format!("fetch {contig}:{ws}-{we}"))?;
        let mut rec = bam::Record::new();
        while let Some(r) = reader.read(&mut rec) {
            r.context("read BAM record")?;
            if rec.flags() & CALLABILITY_FILTER_FLAGS != 0 {
                continue;
            }
            let mq = rec.mapq();
            let (usable, poor) = (mq >= p.min_mapq, mq < p.poor_mapq);
            let mut r0 = rec.pos();
            let mut add = |a: i64, b: i64| {
                let (a, b) = (a.max(ws), b.min(we));
                if a < b {
                    let (i, j) = ((a - ws) as usize, (b - ws) as usize);
                    d_all[i] += 1;
                    d_all[j] -= 1;
                    if usable {
                        d_use[i] += 1;
                        d_use[j] -= 1;
                    }
                    if poor {
                        d_poor[i] += 1;
                        d_poor[j] -= 1;
                    }
                }
            };
            for op in rec.cigar().iter() {
                let l = op.len() as i64;
                match op.char() {
                    'M' | '=' | 'X' | 'D' => {
                        add(r0, r0 + l);
                        r0 += l;
                    }
                    'N' => r0 += l,
                    _ => {}
                }
            }
        }
    }
    let mut hist: Vec<u64> = Vec::new();
    let (mut all, mut usable) = (0i32, 0i32);
    for i in 0..len {
        all += d_all[i];
        usable += d_use[i];
        if all > 0 {
            let u = usable.max(0) as usize;
            if hist.len() <= u {
                hist.resize(u + 1, 0);
            }
            hist[u] += 1;
        }
    }
    let covered: u64 = hist.iter().sum();
    let mut median = 0.0;
    let mut acc = 0u64;
    for (d, c) in hist.iter().enumerate() {
        acc += c;
        if acc * 2 > covered {
            median = d as f64;
            break;
        }
    }
    let excess = p.max_depth_multiple * median;

    let mut out: Vec<CallableInterval> = Vec::new();
    let (mut all, mut usable, mut poor) = (0i32, 0i32, 0i32);
    for i in 0..len {
        all += d_all[i];
        usable += d_use[i];
        poor += d_poor[i];
        let state = if all <= 0 {
            CallState::NoCoverage
        } else if poor as f64 / all as f64 >= p.max_poor_mapq_fraction {
            CallState::PoorMappingQuality
        } else if (usable as u32) < p.min_depth {
            CallState::LowCoverage
        } else if median > 0.0 && usable as f64 > excess {
            CallState::ExcessiveCoverage
        } else {
            CallState::Callable
        };
        let pos = ws as u64 + i as u64;
        match out.last_mut() {
            Some(last) if last.state == state && last.end == pos => last.end = pos + 1,
            _ => out.push(CallableInterval { contig: contig.to_string(), start: pos, end: pos + 1, state }),
        }
    }
    Ok((out, median))
}

pub fn state_at(intervals: &[CallableInterval], contig: &str, pos1: u64) -> Option<CallState> {
    let p = pos1.checked_sub(1)?;
    intervals.iter().find(|iv| iv.contig == contig && iv.start <= p && p < iv.end).map(|iv| iv.state)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn state_lookup() {
        let ivs = vec![
            CallableInterval { contig: "chr1".into(), start: 0, end: 10, state: CallState::LowCoverage },
            CallableInterval { contig: "chr1".into(), start: 10, end: 20, state: CallState::Callable },
        ];
        assert_eq!(state_at(&ivs, "chr1", 1), Some(CallState::LowCoverage));
        assert_eq!(state_at(&ivs, "chr1", 11), Some(CallState::Callable));
        assert_eq!(state_at(&ivs, "chr1", 21), None);
        assert_eq!(state_at(&ivs, "chr2", 5), None);
    }
}
