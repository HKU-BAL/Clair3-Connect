use super::constants::SAMTOOLS_VIEW_FILTER_FLAG;
use rust_htslib::bam::{record::Record, IndexedReader, Read as _};
use rust_htslib::errors::Error as HtsError;
use std::collections::HashSet;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct RawRead {

    pub qname: Vec<u8>,

    pub mq: i8,

    pub raw_mq: u8,

    pub read_start: i64,

    pub read_end: i64,

    pub cigar: Vec<(u32, u32)>,

    pub seqi: Vec<u8>,

    pub qual: Vec<u8>,

    pub l_qseq: u32,

    pub strand: i8,

    pub haplotype: u8,

    pub signal_lengths: Option<Vec<i32>>,
}

pub fn iterate_reads(
    bam_path: &Path,
    region: &(String, i64, i64),
    min_mq: u8,
) -> Result<Vec<RawRead>, HtsError> {
    iterate_reads_inner(bam_path, region, min_mq, false)
}

pub fn iterate_reads_with_dwell(
    bam_path: &Path,
    region: &(String, i64, i64),
    min_mq: u8,
) -> Result<Vec<RawRead>, HtsError> {
    iterate_reads_inner(bam_path, region, min_mq, true)
}

fn iterate_reads_inner(
    bam_path: &Path,
    region: &(String, i64, i64),
    min_mq: u8,
    enable_dwell_time: bool,
) -> Result<Vec<RawRead>, HtsError> {
    let mut out: Vec<RawRead> = Vec::new();
    for_each_read(bam_path, region, min_mq, enable_dwell_time, |r| out.push(r))?;
    Ok(out)
}

pub fn for_each_read<F: FnMut(RawRead)>(
    bam_path: &Path,
    region: &(String, i64, i64),
    min_mq: u8,
    enable_dwell_time: bool,
    mut f: F,
) -> Result<(), HtsError> {
    let mut reader = IndexedReader::from_path(bam_path)?;

    let bgzf_threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
        .min(4);
    if bgzf_threads > 1 {
        let _ = reader.set_threads(bgzf_threads);
    }
    let (contig, start, end) = (region.0.as_str(), region.1, region.2);
    reader.fetch((contig, start, end))?;

    let mut seen_qnames: HashSet<Vec<u8>> = HashSet::new();

    for r in reader.records() {
        let record = r?;
        let core_flag = record.flags();

        if (core_flag & SAMTOOLS_VIEW_FILTER_FLAG) != 0 {
            continue;
        }
        if record.mapq() < min_mq {
            continue;
        }
        let qname = record.qname().to_vec();

        if !seen_qnames.insert(qname.clone()) {
            continue;
        }

        let pos = record.pos();
        let read_end = compute_read_end(pos, &record);
        let cigar = record
            .cigar()
            .iter()
            .map(|c| (c.char() as u32, c.len()))
            .collect();
        let raw_mq = record.mapq();
        let signal_lengths = if enable_dwell_time {
            super::dwell::compute_signal_lengths(&record, record.is_reverse())
        } else {
            None
        };
        f(RawRead {
            qname,
            mq: super::constants::normalize_mq(raw_mq),
            raw_mq,
            read_start: pos,
            read_end,
            cigar,
            seqi: record.seq().encoded.to_vec(),
            qual: record.qual().to_vec(),
            l_qseq: record.seq_len() as u32,
            strand: super::constants::normalize_strand(record.is_reverse()),
            haplotype: super::constants::HAP_UNPHASED,
            signal_lengths,
        });
    }
    Ok(())
}

fn compute_read_end(pos: i64, record: &Record) -> i64 {
    let mut end = pos;
    for c in record.cigar().iter() {

        match c.char() {
            'M' | 'D' | 'N' | '=' | 'X' => end += c.len() as i64,
            _ => {}
        }
    }
    end
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn fixture_bam() -> Option<PathBuf> {
        let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("test_data/ont_chr20_demo/HG003_chr20_demo.bam");
        if p.exists() { Some(p) } else { None }
    }

    #[test]
    fn iterate_filters_smoke() {
        let Some(bam) = fixture_bam() else {
            eprintln!("skipping: no fixture BAM");
            return;
        };

        let region = ("chr20".to_string(), 100_000_i64, 110_000_i64);
        let reads = iterate_reads(&bam, &region, 5).expect("iterate");
        assert!(!reads.is_empty(), "expected at least 1 read in chr20:100k-110k");

        let mut seen = HashSet::new();
        for r in &reads {
            assert!(seen.insert(r.qname.clone()), "duplicate qname leaked");
            assert!(r.mq >= 0, "normalized mq is i8 >= 0");
            assert!(r.read_end >= r.read_start);

            assert!(r.strand == 50 || r.strand == 100);
            assert_eq!(r.haplotype, super::super::constants::HAP_UNPHASED);
        }
    }

    #[test]
    fn filter_flag_drops_secondary() {

        assert_eq!(SAMTOOLS_VIEW_FILTER_FLAG, 2316);

        assert_eq!(0x4 | 0x8 | 0x100 | 0x800, 2316);
    }
}
