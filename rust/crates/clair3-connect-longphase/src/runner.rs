use std::path::{Path, PathBuf};
use std::process::Command;
use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LongPhasePlatform {

    Ont,

    PacBio,
}

impl LongPhasePlatform {
    fn flag(self) -> &'static str {
        match self {
            Self::Ont => "--ont",
            Self::PacBio => "--pb",
        }
    }
}

#[derive(Debug, Error)]
pub enum LongPhaseError {

    #[error("failed to spawn longphase ({binary:?}): {source}")]
    Spawn {

        binary: PathBuf,

        #[source]
        source: std::io::Error,
    },

    #[error("longphase exited with status {status} (stderr: {stderr})")]
    NonZeroExit {

        status: i32,

        stderr: String,
    },

    #[error("longphase succeeded but {path:?} not found")]
    OutputMissing {

        path: PathBuf,
    },
}

#[derive(Debug, Clone)]
pub struct LongPhaseInputs<'a> {

    pub binary: &'a Path,

    pub input_vcf: &'a Path,

    pub bam_path: &'a Path,

    pub fasta_path: &'a Path,

    pub output_prefix: &'a Path,

    pub threads: u32,

    pub platform: LongPhasePlatform,
}

#[derive(Debug, Clone)]
pub struct LongPhaseOutput {

    pub phased_vcf: PathBuf,

    pub stdout: String,

    pub stderr: String,
}

pub fn run_longphase_phase(inputs: &LongPhaseInputs<'_>) -> Result<LongPhaseOutput, LongPhaseError> {
    let mut cmd = Command::new(inputs.binary);
    cmd.arg("phase")
        .arg("-s").arg(inputs.input_vcf)
        .arg("-b").arg(inputs.bam_path)
        .arg("-r").arg(inputs.fasta_path)
        .arg("-t").arg(inputs.threads.to_string())
        .arg("-o").arg(inputs.output_prefix)
        .arg(inputs.platform.flag());

    let output = cmd
        .output()
        .map_err(|e| LongPhaseError::Spawn {
            binary: inputs.binary.to_path_buf(),
            source: e,
        })?;

    let stdout = String::from_utf8_lossy(&output.stdout).into_owned();
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();

    if !output.status.success() {
        return Err(LongPhaseError::NonZeroExit {
            status: output.status.code().unwrap_or(-1),
            stderr,
        });
    }

    let phased_vcf = {
        let mut p = inputs.output_prefix.as_os_str().to_owned();
        p.push(".vcf");
        PathBuf::from(p)
    };
    if !phased_vcf.exists() {
        return Err(LongPhaseError::OutputMissing {
            path: phased_vcf,
        });
    }

    Ok(LongPhaseOutput {
        phased_vcf,
        stdout,
        stderr,
    })
}

pub fn recommended_thread_split(total_threads: u32) -> (u32, u32) {
    let parallel = (total_threads / 2).clamp(1, 4);
    let per_job = (total_threads / parallel).max(1);
    (parallel, per_job)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thread_split_clamps_correctly() {

        assert_eq!(recommended_thread_split(1), (1, 1));

        assert_eq!(recommended_thread_split(4), (2, 2));

        assert_eq!(recommended_thread_split(8), (4, 2));

        assert_eq!(recommended_thread_split(16), (4, 4));

        assert_eq!(recommended_thread_split(32), (4, 8));
    }

    #[test]
    fn platform_flag_strings() {
        assert_eq!(LongPhasePlatform::Ont.flag(), "--ont");
        assert_eq!(LongPhasePlatform::PacBio.flag(), "--pb");
    }

    #[test]
    fn spawn_error_when_binary_missing() {
        let nowhere = PathBuf::from("/nonexistent/__longphase__");
        let r = run_longphase_phase(&LongPhaseInputs {
            binary: &nowhere,
            input_vcf: Path::new("/tmp/x.vcf"),
            bam_path: Path::new("/tmp/x.bam"),
            fasta_path: Path::new("/tmp/x.fa"),
            output_prefix: Path::new("/tmp/x"),
            threads: 1,
            platform: LongPhasePlatform::Ont,
        });
        assert!(matches!(r, Err(LongPhaseError::Spawn { .. })));
    }
}
