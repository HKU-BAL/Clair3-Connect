use anyhow::{anyhow, Context, Result};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::io::Read;
use std::path::{Path, PathBuf};

#[derive(Debug, Deserialize, Serialize)]
pub struct Registry {

    pub server_version: String,

    pub models: Vec<ModelEntry>,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct ModelEntry {

    pub model_id: String,

    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub description: String,

    pub platform: String,

    pub stage: String,

    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub accuracy: Option<String>,

    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub checkpoint_source: Option<CheckpointSource>,

    pub torchscript: TorchScript,

    pub input: TensorSpec,

    pub output: TensorSpec,

    #[serde(default, skip_serializing_if = "serde_yaml::Value::is_null")]
    pub flags: serde_yaml::Value,

    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub notes: Option<String>,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct CheckpointSource {

    pub path: String,

    pub kind: String,

    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sha256: Option<String>,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct TorchScript {

    pub path: String,

    pub sha256: String,
}

#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct TensorSpec {

    pub shape: Vec<i64>,

    pub dtype: String,

    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub layout: Option<String>,
}

impl Registry {

    pub fn load(yaml_path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(yaml_path)
            .with_context(|| format!("read registry {}", yaml_path.display()))?;
        let reg: Self = serde_yaml::from_str(&text).context("parse registry YAML")?;
        Ok(reg)
    }

    pub fn get(&self, model_id: &str) -> Option<&ModelEntry> {
        self.models.iter().find(|m| m.model_id == model_id)
    }

    pub fn verify_torchscript(&self, model: &ModelEntry, registry_dir: &Path) -> Result<PathBuf> {
        let ts_path = registry_dir.join(&model.torchscript.path);
        let actual =
            sha256_file(&ts_path).with_context(|| format!("sha256 of {}", ts_path.display()))?;
        if actual != model.torchscript.sha256 {
            return Err(anyhow!(
                "torchscript sha256 mismatch for {}: file {} has sha256 {}, registry pins {}",
                model.model_id,
                ts_path.display(),
                actual,
                model.torchscript.sha256
            ));
        }
        Ok(ts_path)
    }

    pub fn save(&self, yaml_path: &Path) -> Result<()> {
        let text = serde_yaml::to_string(self).context("serialize registry YAML")?;
        std::fs::write(yaml_path, text)
            .with_context(|| format!("write registry {}", yaml_path.display()))?;
        Ok(())
    }
}

pub fn sha256_file(path: &Path) -> Result<String> {
    let mut hasher = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    let mut fh = std::fs::File::open(path)?;
    loop {
        let n = fh.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(hex::encode(hasher.finalize()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_example_registry() {
        let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();

        let yaml = std::path::Path::new(&manifest_dir)
            .join("..")
            .join("..")
            .join("models")
            .join("model_registry.example.yaml");
        if !yaml.exists() {
            eprintln!("[skip] {} not present", yaml.display());
            return;
        }
        let reg = Registry::load(&yaml).expect("load");
        assert!(!reg.models.is_empty());
        let pileup = reg
            .get("pileup_ont_hac_v500")
            .expect("pileup_ont_hac_v500 not in registry");
        assert_eq!(pileup.stage, "pileup");
        assert_eq!(pileup.input.dtype, "int32");
        assert_eq!(pileup.input.shape, vec![-1, 33, 18]);
        assert_eq!(pileup.output.shape, vec![-1, 24]);
        let fa = reg
            .get("full_alignment_ont_hac_v500")
            .expect("full_alignment_ont_hac_v500 not in registry");
        assert_eq!(fa.stage, "full_alignment");
        assert_eq!(fa.input.dtype, "int8");
        assert_eq!(fa.input.shape, vec![-1, 89, 33, 8]);
    }
}
