#![deny(rust_2018_idioms)]
#![warn(missing_docs)]

pub mod chunk_planner;
pub mod fa_only_pipeline;
pub mod fa_pipeline;
pub mod interactive;
pub mod merge_vcf;
pub mod parallel_driver;
pub mod pileup_pipeline;
pub mod region_filter;
pub mod select_candidates;
pub mod select_het_snp;
pub mod select_qual;
pub mod site_tools;
