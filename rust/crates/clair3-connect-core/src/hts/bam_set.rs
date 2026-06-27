use anyhow::{anyhow, Context, Result};
use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::faidx;
use std::path::{Path, PathBuf};

pub const BAM_PILEUP_FILTER_FLAGS: u16 = 4 + 256 + 2048 + 512 + 1024;

pub struct BamSet {

    pub reader: IndexedReader,

    pub fasta: faidx::Reader,

    pub bam_path: PathBuf,

    pub fasta_path: PathBuf,
}

impl BamSet {

    pub fn open(bam_path: &Path, fasta_path: &Path) -> Result<Self> {
        let reader = IndexedReader::from_path(bam_path)
            .with_context(|| format!("open BAM {}", bam_path.display()))?;
        let fasta = faidx::Reader::from_path(fasta_path)
            .with_context(|| format!("open FASTA {}", fasta_path.display()))?;
        Ok(Self {
            reader,
            fasta,
            bam_path: bam_path.to_path_buf(),
            fasta_path: fasta_path.to_path_buf(),
        })
    }

    pub fn tid(&self, contig: &str) -> Result<u32> {
        let header = self.reader.header();
        let tid = header
            .tid(contig.as_bytes())
            .ok_or_else(|| anyhow!("contig {contig:?} not in BAM header"))?;
        Ok(tid)
    }
}
