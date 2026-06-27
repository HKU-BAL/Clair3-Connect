use std::path::Path;
use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Variant {

    pub position: i64,

    pub ref_base: u8,

    pub alt_base: u8,

    pub genotype: i32,

    pub phase_set: i32,
}

#[derive(Debug, Error)]
pub enum PhasedVcfError {

    #[error("I/O reading {path:?}: {source}")]
    Io {

        path: std::path::PathBuf,

        #[source]
        source: std::io::Error,
    },

    #[error("row {row}: not enough columns (need ≥10, got {got})")]
    MissingColumns {

        row: usize,

        got: usize,
    },
}

pub fn load_phased_vcf(
    path: &Path,
    ctg_filter: Option<&str>,
) -> Result<Vec<Variant>, PhasedVcfError> {
    let text = std::fs::read_to_string(path).map_err(|e| PhasedVcfError::Io {
        path: path.to_path_buf(),
        source: e,
    })?;
    let mut out: Vec<Variant> = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 10 {
            return Err(PhasedVcfError::MissingColumns {
                row: i + 1,
                got: cols.len(),
            });
        }
        let chrom = cols[0];
        if let Some(want) = ctg_filter {
            if chrom != want {
                continue;
            }
        }
        let pos: i64 = match cols[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let ref_field = cols[3];
        let alt_field = cols[4];

        if ref_field.len() != 1 || alt_field.len() != 1 {
            continue;
        }
        let sample = cols[9];
        let format_parts: Vec<&str> = sample.split(':').collect();
        if format_parts.is_empty() {
            continue;
        }
        let genotype_str = format_parts[0];

        if !genotype_str.contains('|') {
            continue;
        }

        let genotype: i32 = if genotype_str == "0|1" { 1 } else { 2 };

        let phase_set: i32 = match format_parts.last().and_then(|s| s.parse().ok()) {
            Some(p) => p,
            None => continue,
        };
        out.push(Variant {
            position: pos - 1,
            ref_base: ref_field.as_bytes()[0],
            alt_base: alt_field.as_bytes()[0],
            genotype,
            phase_set,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_vcf(content: &str) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!(
            "clair3-connect-phased-{}-{}.vcf",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::write(&p, content).unwrap();
        p
    }

    #[test]
    fn parse_simple_phased_snps() {
        let vcf = "##header\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG003\n\
                   chr20\t100\t.\tA\tG\t30\tPASS\t.\tGT:DP:PS\t0|1:30:200\n\
                   chr20\t150\t.\tC\tT\t40\tPASS\t.\tGT:DP:PS\t1|0:25:200\n";
        let p = tmp_vcf(vcf);
        let v = load_phased_vcf(&p, Some("chr20")).unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!(
            v[0],
            Variant {
                position: 99,
                ref_base: b'A',
                alt_base: b'G',
                genotype: 1,
                phase_set: 200,
            }
        );
        assert_eq!(
            v[1],
            Variant {
                position: 149,
                ref_base: b'C',
                alt_base: b'T',
                genotype: 2,
                phase_set: 200,
            }
        );
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn skip_unphased_and_indels() {
        let vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
                   chr1\t10\t.\tA\tG\t30\tPASS\t.\tGT:PS\t0/1:.\n\
                   chr1\t20\t.\tA\tAG\t30\tPASS\t.\tGT:PS\t0|1:200\n\
                   chr1\t30\t.\tA\tG\t30\tPASS\t.\tGT:PS\t1|0:300\n";
        let p = tmp_vcf(vcf);
        let v = load_phased_vcf(&p, None).unwrap();

        assert_eq!(v.len(), 1);
        assert_eq!(v[0].position, 29);
        assert_eq!(v[0].phase_set, 300);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn phased_non_01_genotypes_map_to_2() {

        let vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
                   chr1\t10\t.\tA\tG\t30\tPASS\t.\tGT:PS\t1|1:200\n";
        let p = tmp_vcf(vcf);
        let v = load_phased_vcf(&p, None).unwrap();
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].genotype, 2);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn contig_filter_works() {
        let vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n\
                   chr1\t10\t.\tA\tG\t30\tPASS\t.\tGT:PS\t0|1:200\n\
                   chr2\t20\t.\tA\tG\t30\tPASS\t.\tGT:PS\t0|1:200\n";
        let p = tmp_vcf(vcf);
        let v = load_phased_vcf(&p, Some("chr2")).unwrap();
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].position, 19);
        std::fs::remove_file(p).ok();
    }
}
