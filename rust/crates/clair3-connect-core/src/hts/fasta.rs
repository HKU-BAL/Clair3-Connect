use anyhow::{Context, Result};
use rust_htslib::faidx;

#[derive(Debug, Clone)]
pub struct RefWindow {

    pub seq: Vec<u8>,

    pub start: u64,
}

impl RefWindow {

    pub fn base_at(&self, pos: u64) -> Option<u8> {
        let i = pos.checked_sub(self.start)? as usize;
        self.seq.get(i).copied()
    }
}

pub fn fetch_window(
    fasta: &faidx::Reader,
    contig: &str,
    region_start: u64,
    region_end: u64,
) -> Result<RefWindow> {
    const EXPAND: u64 = 1000;
    let ref_start = region_start.saturating_sub(EXPAND);

    let ref_end_inclusive = region_end.saturating_add(EXPAND).saturating_sub(1);
    let raw = fasta
        .fetch_seq(contig, ref_start as usize, ref_end_inclusive as usize)
        .with_context(|| format!("fetch ref {} [{}-{}]", contig, ref_start, ref_end_inclusive))?;

    let seq: Vec<u8> = raw.iter().map(|b| b.to_ascii_uppercase()).collect();
    Ok(RefWindow {
        seq,
        start: ref_start,
    })
}
