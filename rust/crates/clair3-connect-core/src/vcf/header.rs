use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};
use thiserror::Error;

#[derive(Debug, Clone)]
pub struct VcfHeaderInputs<'a> {

    pub reference_file_path: Option<&'a Path>,

    pub cmdline_file: Option<&'a Path>,

    pub sample_name: &'a str,

    pub clair3_version: &'a str,

    pub gvcf: bool,
}

#[derive(Debug, Clone)]
pub struct VcfHeader {

    pub text: String,

    pub contig_lengths: HashMap<String, u64>,
}

#[derive(Debug, Error)]
pub enum HeaderError {

    #[error("reference index (.fai) not found for FASTA {fa:?}: tried {fai:?}")]
    MissingFai {

        fa: PathBuf,

        fai: PathBuf,
    },

    #[error("failed to read reference index {0:?}: {1}")]
    FaiIo(PathBuf, #[source] io::Error),

    #[error("malformed .fai row {0:?} at line {1}: {2}")]
    FaiMalformed(PathBuf, usize, String),

    #[error("failed to read cmdline file {0:?}: {1}")]
    CmdlineIo(PathBuf, #[source] io::Error),
}

pub fn build_vcf_header(inputs: &VcfHeaderInputs<'_>) -> Result<VcfHeader, HeaderError> {

    let ref_header_str: Option<String> = match inputs.reference_file_path {
        Some(p) if p.exists() => Some(format!("##reference={}", p.display())),
        _ => None,
    };

    let cmdline_str: Option<String> = match inputs.cmdline_file {
        Some(p) if p.exists() => {
            let text = std::fs::read_to_string(p)
                .map_err(|e| HeaderError::CmdlineIo(p.to_path_buf(), e))?;
            let trimmed = text.trim_end();
            if trimmed.is_empty() {
                None
            } else {
                Some(format!("##cmdline={trimmed}"))
            }
        }
        _ => None,
    };

    let base_lines: &[&str] = if inputs.gvcf {
        &[
            "##fileformat=VCFv4.2",
            "##source=Clair3",
            "{clair3_version}",
            r#"##FILTER=<ID=PASS,Description="All filters passed">"#,
            r#"##FILTER=<ID=LowQual,Description="Low quality variant">"#,
            r#"##FILTER=<ID=RefCall,Description="Reference call">"#,
            r#"##INFO=<ID=P,Number=0,Type=Flag,Description="Result from pileup calling">"#,
            r#"##INFO=<ID=F,Number=0,Type=Flag,Description="Result from full-alignment calling">"#,
            r#"##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">"#,
            r#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position (for use with symbolic alleles)">"#,
            r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
            r#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">"#,
            r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads 1. with MQ below 5 or an user-specified threshold, or 2. selected by 'samtools view -F 2316', are filtered)">"#,
            r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">"#,
            r#"##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">"#,
            r#"##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">"#,
            r#"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Observed allele frequency in reads, for each ALT allele, in the same order as listed, or the REF allele for a RefCall">"#,
        ]
    } else {
        &[
            "##fileformat=VCFv4.2",
            "##source=Clair3",
            "{clair3_version}",
            r#"##FILTER=<ID=PASS,Description="All filters passed">"#,
            r#"##FILTER=<ID=LowQual,Description="Low quality variant">"#,
            r#"##FILTER=<ID=RefCall,Description="Reference call">"#,
            r#"##INFO=<ID=P,Number=0,Type=Flag,Description="Result from pileup calling">"#,
            r#"##INFO=<ID=F,Number=0,Type=Flag,Description="Result from full-alignment calling">"#,
            r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#,
            r#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">"#,
            r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads 1. with MQ below 5 or an user-specified threshold, or 2. selected by 'samtools view -F 2316', are filtered)">"#,
            r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">"#,
            r#"##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">"#,
            r#"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Observed allele frequency in reads, for each ALT allele, in the same order as listed, or the REF allele for a RefCall">"#,
        ]
    };

    let mut lines: Vec<String> = base_lines
        .iter()
        .map(|l| {
            if *l == "{clair3_version}" {
                format!("##clair3_version={}", inputs.clair3_version)
            } else {
                (*l).to_string()
            }
        })
        .collect();

    if let Some(ref_line) = &ref_header_str {
        let insert_at = 3.min(lines.len());
        lines.insert(insert_at, ref_line.clone());
    }
    if let Some(cmd_line) = &cmdline_str {
        let insert_at = 3.min(lines.len());
        lines.insert(insert_at, cmd_line.clone());
    }

    let mut text = lines.join("\n");
    text.push('\n');

    let mut contig_lengths: HashMap<String, u64> = HashMap::new();
    if let Some(fa) = inputs.reference_file_path {
        let fai = with_appended_extension(fa, "fai");
        if !fai.exists() {
            return Err(HeaderError::MissingFai {
                fa: fa.to_path_buf(),
                fai,
            });
        }
        let f = File::open(&fai).map_err(|e| HeaderError::FaiIo(fai.clone(), e))?;
        for (i, line) in BufReader::new(f).lines().enumerate() {
            let line = line.map_err(|e| HeaderError::FaiIo(fai.clone(), e))?;
            let mut cols = line.split('\t');
            let name = cols.next().ok_or_else(|| {
                HeaderError::FaiMalformed(fai.clone(), i + 1, line.clone())
            })?;
            let size_str = cols.next().ok_or_else(|| {
                HeaderError::FaiMalformed(fai.clone(), i + 1, line.clone())
            })?;
            let size: u64 = size_str.parse().map_err(|_| {
                HeaderError::FaiMalformed(fai.clone(), i + 1, line.clone())
            })?;
            contig_lengths.insert(name.to_string(), size);
            text.push_str(&format!("##contig=<ID={name},length={size}>\n"));
        }

        text.push_str(&format!(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
            inputs.sample_name
        ));
    }

    Ok(VcfHeader {
        text,
        contig_lengths,
    })
}

fn with_appended_extension(path: &Path, ext: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn unique_tmp(label: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "clair3-connect-{label}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }

    fn make_fixture_fa_fai(contigs: &[(&str, u64)]) -> (PathBuf, PathBuf) {
        let fa = unique_tmp("fa");
        let fai_path = {
            let mut s = fa.as_os_str().to_owned();
            s.push(".fai");
            PathBuf::from(s)
        };
        std::fs::write(&fa, "ACGT").unwrap();
        let mut f = File::create(&fai_path).unwrap();
        for (name, len) in contigs {

            writeln!(f, "{}\t{}\t6\t60\t61", name, len).unwrap();
        }
        (fa, fai_path)
    }

    #[test]
    fn no_reference_emits_only_base_lines_no_chrom() {
        let h = build_vcf_header(&VcfHeaderInputs {
            reference_file_path: None,
            cmdline_file: None,
            sample_name: "HG003",
            clair3_version: "2.0.1",
            gvcf: false,
        })
        .unwrap();
        assert!(h.text.starts_with("##fileformat=VCFv4.2\n"));
        assert!(h.text.contains("##source=Clair3"));
        assert!(h.text.contains("##clair3_version=2.0.1"));
        assert!(
            !h.text.contains("#CHROM"),
            "Python omits #CHROM when reference_file_path=None; got:\n{}",
            h.text
        );
        assert!(h.contig_lengths.is_empty());
    }

    #[test]
    fn version_line_substituted() {
        let h = build_vcf_header(&VcfHeaderInputs {
            reference_file_path: None,
            cmdline_file: None,
            sample_name: "S",
            clair3_version: "9.9.9",
            gvcf: false,
        })
        .unwrap();
        assert!(h.text.contains("##clair3_version=9.9.9"));
        assert!(!h.text.contains("{clair3_version}"));
    }

    #[test]
    fn gvcf_header_extras_present() {
        let h = build_vcf_header(&VcfHeaderInputs {
            reference_file_path: None,
            cmdline_file: None,
            sample_name: "S",
            clair3_version: "2.0.1",
            gvcf: true,
        })
        .unwrap();
        assert!(h.text.contains("##ALT=<ID=NON_REF"));
        assert!(h.text.contains("##INFO=<ID=END"));
        assert!(h.text.contains("MIN_DP"));
    }

    #[test]
    fn reference_emits_chrom_and_contig_lines() {
        let (fa, fai) = make_fixture_fa_fai(&[("chr20", 64444167)]);
        let h = build_vcf_header(&VcfHeaderInputs {
            reference_file_path: Some(&fa),
            cmdline_file: None,
            sample_name: "S",
            clair3_version: "2.0.1",
            gvcf: false,
        })
        .unwrap();
        assert!(h.text.contains(&format!("##reference={}", fa.display())));
        assert!(h.text.contains("##contig=<ID=chr20,length=64444167>"));
        assert!(h
            .text
            .contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"));
        assert_eq!(h.contig_lengths.get("chr20"), Some(&64444167));
        std::fs::remove_file(fa).ok();
        std::fs::remove_file(fai).ok();
    }

    #[test]
    fn missing_fai_is_an_error() {
        let fa = unique_tmp("fa");
        std::fs::write(&fa, "ACGT").unwrap();
        let r = build_vcf_header(&VcfHeaderInputs {
            reference_file_path: Some(&fa),
            cmdline_file: None,
            sample_name: "S",
            clair3_version: "2.0.1",
            gvcf: false,
        });
        assert!(matches!(r, Err(HeaderError::MissingFai { .. })));
        std::fs::remove_file(fa).ok();
    }

    #[test]
    fn nonexistent_reference_path_triggers_missing_fai() {

        let fake = PathBuf::from("/nonexistent/__missing__.fa");
        let r = build_vcf_header(&VcfHeaderInputs {
            reference_file_path: Some(&fake),
            cmdline_file: None,
            sample_name: "S",
            clair3_version: "2.0.1",
            gvcf: false,
        });
        assert!(
            matches!(r, Err(HeaderError::MissingFai { .. })),
            "expected MissingFai, got {r:?}"
        );
    }

    #[test]
    fn cmdline_file_contents_used() {
        let cmd_file = unique_tmp("cmd");
        std::fs::write(&cmd_file, "clair3 --bam x.bam --ref y.fa\n   \n").unwrap();
        let (fa, fai) = make_fixture_fa_fai(&[("chr1", 1000)]);
        let h = build_vcf_header(&VcfHeaderInputs {
            reference_file_path: Some(&fa),
            cmdline_file: Some(&cmd_file),
            sample_name: "S",
            clair3_version: "2.0.1",
            gvcf: false,
        })
        .unwrap();

        assert!(
            h.text.contains("##cmdline=clair3 --bam x.bam --ref y.fa"),
            "got:\n{}",
            h.text
        );
        std::fs::remove_file(fa).ok();
        std::fs::remove_file(fai).ok();
        std::fs::remove_file(cmd_file).ok();
    }

    #[test]
    fn cmdline_file_missing_or_empty_omits_line() {
        let missing = PathBuf::from("/nonexistent/__cmd__.txt");
        let (fa, fai) = make_fixture_fa_fai(&[("chr1", 100)]);
        let h = build_vcf_header(&VcfHeaderInputs {
            reference_file_path: Some(&fa),
            cmdline_file: Some(&missing),
            sample_name: "S",
            clair3_version: "2.0.1",
            gvcf: false,
        })
        .unwrap();
        assert!(!h.text.contains("##cmdline="));

        let empty = unique_tmp("emptycmd");
        std::fs::write(&empty, "   \n").unwrap();
        let h2 = build_vcf_header(&VcfHeaderInputs {
            reference_file_path: Some(&fa),
            cmdline_file: Some(&empty),
            sample_name: "S",
            clair3_version: "2.0.1",
            gvcf: false,
        })
        .unwrap();
        assert!(!h2.text.contains("##cmdline="));

        std::fs::remove_file(fa).ok();
        std::fs::remove_file(fai).ok();
        std::fs::remove_file(empty).ok();
    }

    #[test]
    fn cmdline_and_reference_both_at_index_3() {
        let cmd_file = unique_tmp("cmd2");
        std::fs::write(&cmd_file, "the cmd").unwrap();
        let (fa, fai) = make_fixture_fa_fai(&[("chr1", 100)]);
        let h = build_vcf_header(&VcfHeaderInputs {
            reference_file_path: Some(&fa),
            cmdline_file: Some(&cmd_file),
            sample_name: "S",
            clair3_version: "2.0.1",
            gvcf: false,
        })
        .unwrap();
        let lines: Vec<&str> = h.text.lines().collect();
        assert_eq!(lines[3], "##cmdline=the cmd");
        assert_eq!(lines[4], format!("##reference={}", fa.display()));
        std::fs::remove_file(fa).ok();
        std::fs::remove_file(fai).ok();
        std::fs::remove_file(cmd_file).ok();
    }

    #[test]
    fn returns_contig_lengths_map() {
        let (fa, fai) = make_fixture_fa_fai(&[("chr1", 100), ("chr2", 200), ("chrX", 999)]);
        let h = build_vcf_header(&VcfHeaderInputs {
            reference_file_path: Some(&fa),
            cmdline_file: None,
            sample_name: "S",
            clair3_version: "2.0.1",
            gvcf: false,
        })
        .unwrap();
        assert_eq!(h.contig_lengths.get("chr1"), Some(&100));
        assert_eq!(h.contig_lengths.get("chr2"), Some(&200));
        assert_eq!(h.contig_lengths.get("chrX"), Some(&999));
        std::fs::remove_file(fa).ok();
        std::fs::remove_file(fai).ok();
    }
}
