use serde::{Deserialize, Serialize};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use thiserror::Error;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ManifestEntry {

    pub chunk_id: String,

    pub contig: String,

    pub start_1based: u64,

    pub end_1based: u64,

    pub row_count: u64,

    pub vcf_sha256: String,

    pub vcf_path: String,

    pub source_phase: String,

    pub completed_at_utc: String,
}

impl ManifestEntry {

    pub fn make_chunk_id(contig: &str, start_1based: u64, end_1based: u64) -> String {
        format!("{}:{}-{}", contig, start_1based, end_1based)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ManifestWarning {

    MalformedLine {

        line_number: usize,

        snippet: String,
    },
}

#[derive(Debug, Error)]
pub enum ManifestError {

    #[error("manifest I/O: {0}")]
    Io(#[from] std::io::Error),

    #[error("manifest JSONL decode: {0}")]
    Json(#[from] serde_json::Error),

    #[error("chunk {chunk_id:?} VCF sha256 mismatch (expected {expected}, got {actual})")]
    ShaMismatch {

        chunk_id: String,

        expected: String,

        actual: String,
    },
}

pub struct ManifestWriter {
    dir: PathBuf,
    file: Mutex<BufWriter<File>>,
}

impl ManifestWriter {

    pub fn open(dir: &Path) -> Result<Self, ManifestError> {
        std::fs::create_dir_all(dir)?;
        let path = dir.join("manifest.jsonl");
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)?;
        Ok(Self {
            dir: dir.to_path_buf(),
            file: Mutex::new(BufWriter::new(file)),
        })
    }

    pub fn record_chunk(&self, entry: &ManifestEntry) -> Result<(), ManifestError> {
        let json = serde_json::to_string(entry)?;
        let mut guard = self.file.lock().expect("manifest mutex");
        guard.write_all(json.as_bytes())?;
        guard.write_all(b"\n")?;
        guard.flush()?;
        guard.get_ref().sync_data()?;
        Ok(())
    }

    pub fn store_chunk_vcf(&self, vcf_path_rel: &str, vcf_text: &str) -> Result<String, ManifestError> {
        let full = self.dir.join(vcf_path_rel);
        if let Some(parent) = full.parent() {
            std::fs::create_dir_all(parent)?;
        }
        let mut f = OpenOptions::new()
            .create(true)
            .truncate(true)
            .write(true)
            .open(&full)?;
        f.write_all(vcf_text.as_bytes())?;
        f.flush()?;
        f.sync_all()?;
        Ok(sha256_hex(vcf_text.as_bytes()))
    }

    pub fn dir(&self) -> &Path {
        &self.dir
    }
}

pub fn read_manifest(dir: &Path) -> Result<Vec<ManifestEntry>, ManifestError> {
    Ok(read_manifest_with_warnings(dir)?.0)
}

pub fn read_manifest_with_warnings(
    dir: &Path,
) -> Result<(Vec<ManifestEntry>, Vec<ManifestWarning>), ManifestError> {
    let path = dir.join("manifest.jsonl");
    if !path.exists() {
        return Ok((Vec::new(), Vec::new()));
    }
    let f = File::open(&path)?;
    let r = BufReader::new(f);
    let mut entries = Vec::new();
    let mut warnings = Vec::new();
    for (i, line) in r.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        match serde_json::from_str::<ManifestEntry>(&line) {
            Ok(e) => entries.push(e),
            Err(_) => {
                let snippet: String = line.chars().take(200).collect();
                warnings.push(ManifestWarning::MalformedLine {
                    line_number: i + 1,
                    snippet,
                });
                break;
            }
        }
    }
    Ok((entries, warnings))
}

pub fn verify_manifest(dir: &Path, entries: &[ManifestEntry]) -> Result<(), ManifestError> {
    for e in entries {
        let path = dir.join(&e.vcf_path);
        if !path.exists() {
            return Err(ManifestError::Io(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!("chunk vcf missing: {}", path.display()),
            )));
        }
        let bytes = std::fs::read(&path)?;
        let actual = sha256_hex(&bytes);
        if actual != e.vcf_sha256 {
            return Err(ManifestError::ShaMismatch {
                chunk_id: e.chunk_id.clone(),
                expected: e.vcf_sha256.clone(),
                actual,
            });
        }
    }
    Ok(())
}

fn sha256_hex(bytes: &[u8]) -> String {
    use sha2::{Digest, Sha256};
    let mut h = Sha256::new();
    h.update(bytes);
    let v = h.finalize();
    let mut s = String::with_capacity(64);
    for b in v.iter() {
        s.push_str(&format!("{:02x}", b));
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_dir(label: &str) -> PathBuf {
        let p = std::env::temp_dir().join(format!(
            "clair3-connect-manifest-{label}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&p).unwrap();
        p
    }

    #[test]
    fn write_and_read_roundtrip() {
        let dir = tmp_dir("rt");
        let w = ManifestWriter::open(&dir).unwrap();
        let vcf_text = "chr20\t100\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n";
        let sha = w.store_chunk_vcf("chunk-00000.vcf", vcf_text).unwrap();
        let entry = ManifestEntry {
            chunk_id: ManifestEntry::make_chunk_id("chr20", 0, 1),
            contig: "chr20".to_string(),
            start_1based: 1,
            end_1based: 1_000_000,
            row_count: 1,
            vcf_sha256: sha,
            vcf_path: "chunk-00000.vcf".to_string(),
            source_phase: "pileup".to_string(),
            completed_at_utc: "2026-05-14T00:00:00Z".to_string(),
        };
        w.record_chunk(&entry).unwrap();
        drop(w);

        let entries = read_manifest(&dir).unwrap();
        assert_eq!(entries, vec![entry]);
        verify_manifest(&dir, &entries).unwrap();
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn partial_last_line_is_dropped() {
        let dir = tmp_dir("partial");
        let w = ManifestWriter::open(&dir).unwrap();
        let vcf_text = "row\n";
        let sha = w.store_chunk_vcf("chunk-00000.vcf", vcf_text).unwrap();
        w.record_chunk(&ManifestEntry {
            chunk_id: "chr1:0/1".to_string(),
            contig: "chr1".to_string(),
            start_1based: 1,
            end_1based: 100,
            row_count: 1,
            vcf_sha256: sha,
            vcf_path: "chunk-00000.vcf".to_string(),
            source_phase: "pileup".to_string(),
            completed_at_utc: "2026-05-14T00:00:00Z".to_string(),
        })
        .unwrap();
        drop(w);

        let path = dir.join("manifest.jsonl");
        let mut f = OpenOptions::new().append(true).open(&path).unwrap();
        f.write_all(b"{\"chunk_id\":\"chr1:1/").unwrap();
        drop(f);

        let entries = read_manifest(&dir).unwrap();
        assert_eq!(entries.len(), 1, "partial last line should be ignored");
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn verify_catches_sha_mismatch() {
        let dir = tmp_dir("sha");
        let w = ManifestWriter::open(&dir).unwrap();
        let vcf_text = "row\n";
        let sha = w.store_chunk_vcf("chunk-00000.vcf", vcf_text).unwrap();
        w.record_chunk(&ManifestEntry {
            chunk_id: "chr1:0/1".to_string(),
            contig: "chr1".to_string(),
            start_1based: 1,
            end_1based: 100,
            row_count: 1,
            vcf_sha256: sha,
            vcf_path: "chunk-00000.vcf".to_string(),
            source_phase: "pileup".to_string(),
            completed_at_utc: "2026-05-14T00:00:00Z".to_string(),
        })
        .unwrap();
        drop(w);

        std::fs::write(dir.join("chunk-00000.vcf"), "DIFFERENT").unwrap();
        let entries = read_manifest(&dir).unwrap();
        let r = verify_manifest(&dir, &entries);
        assert!(matches!(r, Err(ManifestError::ShaMismatch { .. })));
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn missing_manifest_returns_empty() {
        let dir = tmp_dir("missing");

        let entries = read_manifest(&dir).unwrap();
        assert!(entries.is_empty());
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn chunk_id_format_is_stable() {

        assert_eq!(
            ManifestEntry::make_chunk_id("chr20", 1, 5_000_000),
            "chr20:1-5000000"
        );
        assert_eq!(
            ManifestEntry::make_chunk_id("chr20", 5_000_001, 10_000_000),
            "chr20:5000001-10000000"
        );
    }

    #[test]
    fn read_with_warnings_surfaces_partial_line() {
        let dir = tmp_dir("warn");
        let w = ManifestWriter::open(&dir).unwrap();
        let vcf = "row\n";
        let sha = w.store_chunk_vcf("chunk-00000.vcf", vcf).unwrap();
        w.record_chunk(&ManifestEntry {
            chunk_id: ManifestEntry::make_chunk_id("chr1", 1, 100),
            contig: "chr1".to_string(),
            start_1based: 1,
            end_1based: 100,
            row_count: 1,
            vcf_sha256: sha,
            vcf_path: "chunk-00000.vcf".to_string(),
            source_phase: "pileup".to_string(),
            completed_at_utc: "2026-05-14T00:00:00Z".to_string(),
        })
        .unwrap();
        drop(w);

        let path = dir.join("manifest.jsonl");
        let mut f = OpenOptions::new().append(true).open(&path).unwrap();
        f.write_all(b"{\"chunk_id\":\"bad").unwrap();
        drop(f);

        let (entries, warnings) = read_manifest_with_warnings(&dir).unwrap();
        assert_eq!(entries.len(), 1);
        assert_eq!(warnings.len(), 1);
        match &warnings[0] {
            ManifestWarning::MalformedLine { line_number, snippet } => {
                assert_eq!(*line_number, 2);
                assert!(snippet.contains("bad"));
            }
        }
        std::fs::remove_dir_all(dir).ok();
    }
}
