#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChunkRegion {

    pub contig: String,

    pub start_1based: u64,

    pub end_1based: u64,

    pub chunk_idx: u32,

    pub chunk_count: u32,
}

pub const DEFAULT_CHUNK_SIZE: u64 = 5_000_000;

pub const MIN_CHUNK_SIZE: u64 = 1_000;

pub fn plan_chunks(
    contig: &str,
    contig_length: u64,
    chunk_size_bp: u64,
) -> Vec<ChunkRegion> {
    if contig_length == 0 {
        return Vec::new();
    }
    assert!(
        chunk_size_bp >= MIN_CHUNK_SIZE,
        "chunk_size_bp ({}) must be >= {} to avoid pathological chunk counts",
        chunk_size_bp,
        MIN_CHUNK_SIZE
    );
    let chunk_size = chunk_size_bp;
    let chunk_count_full = contig_length / chunk_size;
    let remainder = contig_length % chunk_size;
    let chunk_count = if remainder > 0 {
        chunk_count_full + 1
    } else {
        chunk_count_full
    }
    .max(1);

    let mut out = Vec::with_capacity(chunk_count as usize);
    for i in 0..chunk_count {
        let start = i * chunk_size + 1;
        let end = ((i + 1) * chunk_size).min(contig_length);
        out.push(ChunkRegion {
            contig: contig.to_string(),
            start_1based: start,
            end_1based: end,
            chunk_idx: i as u32,
            chunk_count: chunk_count as u32,
        });
    }
    out
}

pub fn plan_chunks_from_fai(
    fai_entries: &[(String, u64)],
    chunk_size_bp: u64,
) -> Vec<ChunkRegion> {
    let mut out = Vec::new();
    for (name, length) in fai_entries {
        out.extend(plan_chunks(name, *length, chunk_size_bp));
    }
    out
}

pub fn read_fai(fai_path: &std::path::Path) -> std::io::Result<Vec<(String, u64)>> {
    let text = std::fs::read_to_string(fai_path)?;
    let mut out = Vec::new();
    for line in text.lines() {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 {
            continue;
        }
        let name = cols[0].to_string();
        let length: u64 = match cols[1].parse() {
            Ok(n) => n,
            Err(_) => continue,
        };
        out.push((name, length));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_contig_gets_one_chunk() {
        let r = plan_chunks("chr1", 1_000, 5_000_000);
        assert_eq!(r.len(), 1);
        assert_eq!(r[0].start_1based, 1);
        assert_eq!(r[0].end_1based, 1_000);
        assert_eq!(r[0].chunk_idx, 0);
        assert_eq!(r[0].chunk_count, 1);
    }

    #[test]
    fn even_multiple_of_chunk_size() {
        let r = plan_chunks("chr1", 15_000_000, 5_000_000);
        assert_eq!(r.len(), 3);
        assert_eq!(r[0].start_1based, 1);
        assert_eq!(r[0].end_1based, 5_000_000);
        assert_eq!(r[1].start_1based, 5_000_001);
        assert_eq!(r[1].end_1based, 10_000_000);
        assert_eq!(r[2].start_1based, 10_000_001);
        assert_eq!(r[2].end_1based, 15_000_000);
    }

    #[test]
    fn remainder_produces_extra_chunk() {
        let r = plan_chunks("chr1", 12_345_678, 5_000_000);
        assert_eq!(r.len(), 3);
        assert_eq!(r[2].start_1based, 10_000_001);
        assert_eq!(r[2].end_1based, 12_345_678);
    }

    #[test]
    fn zero_length_contig_returns_empty() {
        let r = plan_chunks("chr1", 0, 5_000_000);
        assert!(r.is_empty());
    }

    #[test]
    fn chr20_demo_one_chunk() {

        let r = plan_chunks("chr20", 64_444_167, DEFAULT_CHUNK_SIZE);
        assert_eq!(r.len(), 13);
        assert_eq!(r[12].end_1based, 64_444_167);
    }

    #[test]
    #[should_panic(expected = "chunk_size_bp")]
    fn plan_chunks_rejects_tiny_chunk_size() {
        let _ = plan_chunks("chr1", 1_000_000, 0);
    }

    #[test]
    #[should_panic(expected = "chunk_size_bp")]
    fn plan_chunks_rejects_below_min() {
        let _ = plan_chunks("chr1", 1_000_000, 100);
    }

    #[test]
    fn read_fai_basic() {
        let tmp = std::env::temp_dir().join(format!(
            "clair3-connect-fai-{}.fai",
            std::process::id()
        ));
        std::fs::write(&tmp, "chr1\t1000\t6\t60\t61\nchr2\t2000\t1100\t60\t61\n").unwrap();
        let entries = read_fai(&tmp).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0], ("chr1".to_string(), 1000));
        assert_eq!(entries[1], ("chr2".to_string(), 2000));
        std::fs::remove_file(tmp).ok();
    }
}
