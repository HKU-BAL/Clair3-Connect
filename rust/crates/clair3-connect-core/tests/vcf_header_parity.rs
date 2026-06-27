use clair3_connect_core::vcf::header::{build_vcf_header, VcfHeaderInputs};
use std::path::PathBuf;

const FIXTURE_REF_FA: &str =
    "<PROJECT_ROOT>/test_data/ont_chr20_demo/GRCh38_no_alt_chr20.fa";

const EXPECTED_WITH_PRINT_NEWLINE: &str = concat!(
    "##fileformat=VCFv4.2\n",
    "##source=Clair3\n",
    "##clair3_version=2.0.1\n",
    "##reference=<PROJECT_ROOT>/test_data/ont_chr20_demo/GRCh38_no_alt_chr20.fa\n",
    "##FILTER=<ID=PASS,Description=\"All filters passed\">\n",
    "##FILTER=<ID=LowQual,Description=\"Low quality variant\">\n",
    "##FILTER=<ID=RefCall,Description=\"Reference call\">\n",
    "##INFO=<ID=P,Number=0,Type=Flag,Description=\"Result from pileup calling\">\n",
    "##INFO=<ID=F,Number=0,Type=Flag,Description=\"Result from full-alignment calling\">\n",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads 1. with MQ below 5 or an user-specified threshold, or 2. selected by 'samtools view -F 2316', are filtered)\">\n",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n",
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n",
    "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Observed allele frequency in reads, for each ALT allele, in the same order as listed, or the REF allele for a RefCall\">\n",
    "##contig=<ID=chr20,length=64444167>\n",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG003\n",
);

#[test]
fn header_byte_identical_to_python_for_fixture() {
    let fa = PathBuf::from(FIXTURE_REF_FA);
    if !fa.exists() {
        eprintln!("Skipping — fixture FASTA missing at {}", FIXTURE_REF_FA);
        return;
    }
    let got = build_vcf_header(&VcfHeaderInputs {
        reference_file_path: Some(&fa),
        cmdline_file: None,
        sample_name: "HG003",
        clair3_version: "2.0.1",
        gvcf: false,
    })
    .expect("build_vcf_header");

    let expected_raw = EXPECTED_WITH_PRINT_NEWLINE.strip_suffix('\n').unwrap();

    if got.text != expected_raw {
        let g_lines: Vec<&str> = got.text.lines().collect();
        let e_lines: Vec<&str> = expected_raw.lines().collect();
        for (i, (g, e)) in g_lines.iter().zip(e_lines.iter()).enumerate() {
            if g != e {
                panic!("line {i} differs:\n got: {g}\n exp: {e}");
            }
        }
        panic!(
            "line counts differ: got={} exp={}",
            g_lines.len(),
            e_lines.len()
        );
    }

    assert_eq!(got.contig_lengths.get("chr20"), Some(&64444167));
}
