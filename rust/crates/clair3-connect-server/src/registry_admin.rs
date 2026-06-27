use anyhow::{bail, Context, Result};
use clair3_connect_infer::model_registry::{
    sha256_file, ModelEntry, Registry, TensorSpec, TorchScript,
};
use console::style;
use dialoguer::{theme::ColorfulTheme, Confirm, Input, Select};
use std::path::{Path, PathBuf};

const PLATFORMS: &[&str] = &["ont", "hifi", "ilmn"];
const STAGES: &[&str] = &["pileup", "full_alignment"];

pub fn print_registry_status(yaml: &Path) {
    print_registry_summary(&inspect_registry(yaml));
}

pub fn manage(theme: &ColorfulTheme, yaml: &Path) -> Result<()> {
    let registry_dir = yaml
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."));

    let mut reg = match Registry::load(yaml) {
        Ok(r) => r,
        Err(e) => {
            println!("  {}: {e:#}", style("cannot load registry").red());
            return Ok(());
        }
    };
    let mut dirty = false;

    loop {
        let items = [
            "List entries (verbose)",
            "Verify all entries (sha + file presence)",
            "Add a model",
            "Edit a model",
            "Remove a model",
            "Save changes to YAML",
            if dirty { "Done (unsaved changes!)" } else { "Done" },
        ];
        let choice = Select::with_theme(theme)
            .with_prompt(format!(
                "  Model registry — {} entr{} loaded{}",
                reg.models.len(),
                if reg.models.len() == 1 { "y" } else { "ies" },
                if dirty { " (modified)" } else { "" }
            ))
            .items(&items)
            .default(if dirty { 5 } else { 6 })
            .interact()
            .context("registry menu")?;
        match choice {
            0 => print_entries_verbose(&reg),
            1 => print_registry_summary(&inspect_registry(yaml)),
            2 => {
                let entry = prompt_model_entry(theme, &registry_dir, None)?;
                if reg.models.iter().any(|m| m.model_id == entry.model_id) {
                    println!(
                        "  {}: a model with id '{}' already exists — not added",
                        style("error").red(),
                        entry.model_id
                    );
                } else {
                    println!("  {} added '{}'", style("✓").green(), entry.model_id);
                    reg.models.push(entry);
                    dirty = true;
                }
            }
            3 => {
                if let Some(i) = pick_entry(theme, &reg, "Edit which model?")? {
                    let edited =
                        prompt_model_entry(theme, &registry_dir, Some(&reg.models[i]))?;

                    if reg
                        .models
                        .iter()
                        .enumerate()
                        .any(|(j, m)| j != i && m.model_id == edited.model_id)
                    {
                        println!(
                            "  {}: id '{}' collides with another entry — edit discarded",
                            style("error").red(),
                            edited.model_id
                        );
                    } else {
                        println!("  {} updated '{}'", style("✓").green(), edited.model_id);
                        reg.models[i] = edited;
                        dirty = true;
                    }
                }
            }
            4 => {
                if let Some(i) = pick_entry(theme, &reg, "Remove which model?")? {
                    let id = reg.models[i].model_id.clone();
                    let yes = Confirm::with_theme(theme)
                        .with_prompt(format!("  Really remove '{id}'?"))
                        .default(false)
                        .interact()
                        .context("confirm remove")?;
                    if yes {
                        reg.models.remove(i);
                        dirty = true;
                        println!("  {} removed '{}'", style("✓").green(), id);
                    }
                }
            }
            5 => {
                save_registry(&reg, yaml)?;
                dirty = false;
            }
            6 => {
                if dirty {
                    let save = Confirm::with_theme(theme)
                        .with_prompt("  You have unsaved changes — save before leaving?")
                        .default(true)
                        .interact()
                        .context("confirm save-on-exit")?;
                    if save {
                        save_registry(&reg, yaml)?;
                    }
                }
                break;
            }
            _ => unreachable!(),
        }
    }
    Ok(())
}

fn save_registry(reg: &Registry, yaml: &Path) -> Result<()> {
    if yaml.exists() {
        let bak = yaml.with_extension("yaml.bak");
        std::fs::copy(yaml, &bak)
            .with_context(|| format!("back up {} -> {}", yaml.display(), bak.display()))?;
        println!("  {} backed up original to {}", style("·").dim(), bak.display());
    }
    reg.save(yaml).context("save registry")?;
    println!(
        "  {} wrote {} ({} entries)",
        style("✓").green(),
        yaml.display(),
        reg.models.len()
    );
    Ok(())
}

fn pick_entry(
    theme: &ColorfulTheme,
    reg: &Registry,
    prompt: &str,
) -> Result<Option<usize>> {
    if reg.models.is_empty() {
        println!("  {}", style("(registry has no entries)").yellow());
        return Ok(None);
    }
    let mut labels: Vec<String> = reg
        .models
        .iter()
        .map(|m| format!("{} [{} / {}]", m.model_id, m.platform, m.stage))
        .collect();
    labels.push("(cancel)".to_string());
    let idx = Select::with_theme(theme)
        .with_prompt(format!("  {prompt}"))
        .items(&labels)
        .default(labels.len() - 1)
        .interact()
        .context("pick entry")?;
    Ok(if idx + 1 == labels.len() {
        None
    } else {
        Some(idx)
    })
}

fn prompt_model_entry(
    theme: &ColorfulTheme,
    registry_dir: &Path,
    existing: Option<&ModelEntry>,
) -> Result<ModelEntry> {
    println!();
    println!(
        "  {}",
        style(if existing.is_some() {
            "Edit model — press Enter to keep each current value"
        } else {
            "Add model — fill in each field"
        })
        .bold()
    );

    let model_id: String = Input::with_theme(theme)
        .with_prompt("  model_id")
        .with_initial_text(existing.map(|e| e.model_id.clone()).unwrap_or_default())
        .interact_text()
        .context("read model_id")?;

    let platform = pick_or_default(theme, "  platform", PLATFORMS, existing.map(|e| e.platform.as_str()))?;
    let stage = pick_or_default(theme, "  stage", STAGES, existing.map(|e| e.stage.as_str()))?;

    let accuracy: String = Input::with_theme(theme)
        .with_prompt("  accuracy (informational, blank = none)")
        .allow_empty(true)
        .with_initial_text(
            existing
                .and_then(|e| e.accuracy.clone())
                .unwrap_or_default(),
        )
        .interact_text()
        .context("read accuracy")?;

    let description: String = Input::with_theme(theme)
        .with_prompt("  description (blank = none)")
        .allow_empty(true)
        .with_initial_text(existing.map(|e| e.description.clone()).unwrap_or_default())
        .interact_text()
        .context("read description")?;

    let ts_path: String = Input::with_theme(theme)
        .with_prompt("  torchscript path (relative to registry dir)")
        .with_initial_text(
            existing
                .map(|e| e.torchscript.path.clone())
                .unwrap_or_default(),
        )
        .interact_text()
        .context("read torchscript path")?;
    let ts_sha = resolve_ts_sha(
        theme,
        registry_dir,
        &ts_path,
        existing.map(|e| e.torchscript.sha256.as_str()),
    )?;

    let input = prompt_tensor_spec(theme, "input", existing.map(|e| &e.input))?;
    let output = prompt_tensor_spec(theme, "output", existing.map(|e| &e.output))?;
    let flags = prompt_flags(theme, existing.map(|e| &e.flags))?;

    Ok(ModelEntry {
        model_id,
        description,
        platform,
        stage,
        accuracy: if accuracy.is_empty() { None } else { Some(accuracy) },

        checkpoint_source: existing.and_then(|e| e.checkpoint_source.clone()),
        torchscript: TorchScript { path: ts_path, sha256: ts_sha },
        input,
        output,
        flags,
        notes: existing.and_then(|e| e.notes.clone()),
    })
}

fn resolve_ts_sha(
    theme: &ColorfulTheme,
    registry_dir: &Path,
    ts_path: &str,
    previous: Option<&str>,
) -> Result<String> {
    let resolved = registry_dir.join(ts_path);
    if resolved.exists() {
        match sha256_file(&resolved) {
            Ok(sha) => {
                println!(
                    "  {} computed sha256 of {} = {}…",
                    style("·").dim(),
                    resolved.display(),
                    &sha[..sha.len().min(16)]
                );
                return Ok(sha);
            }
            Err(e) => println!("  {}: {e:#}", style("sha256 failed").yellow()),
        }
    } else {
        println!(
            "  {}: {} not found — the server will treat this as a template entry",
            style("warning").yellow(),
            resolved.display()
        );
    }
    let typed: String = Input::with_theme(theme)
        .with_prompt("  torchscript sha256 (64 hex; blank = all-zero placeholder)")
        .allow_empty(true)
        .with_initial_text(previous.unwrap_or("").to_string())
        .interact_text()
        .context("read sha256")?;
    Ok(if typed.is_empty() {
        "0".repeat(64)
    } else {
        typed
    })
}

fn pick_or_default(
    theme: &ColorfulTheme,
    prompt: &str,
    options: &[&str],
    current: Option<&str>,
) -> Result<String> {
    let default_idx = current
        .and_then(|c| options.iter().position(|o| *o == c))
        .unwrap_or(0);
    let idx = Select::with_theme(theme)
        .with_prompt(prompt)
        .items(options)
        .default(default_idx)
        .interact()
        .with_context(|| format!("select {prompt}"))?;
    Ok(options[idx].to_string())
}

fn prompt_tensor_spec(
    theme: &ColorfulTheme,
    which: &str,
    existing: Option<&TensorSpec>,
) -> Result<TensorSpec> {
    let shape_default = existing
        .map(|t| {
            t.shape
                .iter()
                .map(|d| d.to_string())
                .collect::<Vec<_>>()
                .join(",")
        })
        .unwrap_or_default();
    let shape_str: String = Input::with_theme(theme)
        .with_prompt(format!("  {which} shape (comma-separated, -1 = dynamic batch)"))
        .with_initial_text(shape_default)
        .interact_text()
        .with_context(|| format!("read {which} shape"))?;
    let shape = parse_shape(&shape_str)
        .with_context(|| format!("parse {which} shape {shape_str:?}"))?;

    let dtype: String = Input::with_theme(theme)
        .with_prompt(format!("  {which} dtype (int8 | int32 | float32 | …)"))
        .with_initial_text(
            existing
                .map(|t| t.dtype.clone())
                .unwrap_or_else(|| "int32".to_string()),
        )
        .interact_text()
        .with_context(|| format!("read {which} dtype"))?;

    let layout: String = Input::with_theme(theme)
        .with_prompt(format!("  {which} layout hint (blank = none, e.g. NHWC)"))
        .allow_empty(true)
        .with_initial_text(
            existing
                .and_then(|t| t.layout.clone())
                .unwrap_or_default(),
        )
        .interact_text()
        .with_context(|| format!("read {which} layout"))?;

    Ok(TensorSpec {
        shape,
        dtype,
        layout: if layout.is_empty() { None } else { Some(layout) },
    })
}

fn parse_shape(s: &str) -> Result<Vec<i64>> {
    let mut out = Vec::new();
    for tok in s.split(',') {
        let t = tok.trim();
        if t.is_empty() {
            continue;
        }
        out.push(t.parse::<i64>().with_context(|| format!("bad dim {t:?}"))?);
    }
    if out.is_empty() {
        bail!("shape must have at least one dimension");
    }
    Ok(out)
}

fn prompt_flags(
    theme: &ColorfulTheme,
    existing: Option<&serde_yaml::Value>,
) -> Result<serde_yaml::Value> {
    let get_bool = |key: &str, default: bool| -> bool {
        existing
            .and_then(|v| v.get(key))
            .and_then(serde_yaml::Value::as_bool)
            .unwrap_or(default)
    };
    let add_indel = Confirm::with_theme(theme)
        .with_prompt("  flag add_indel_length")
        .default(get_bool("add_indel_length", false))
        .interact()
        .context("flag add_indel_length")?;
    let predict = Confirm::with_theme(theme)
        .with_prompt("  flag predict")
        .default(get_bool("predict", true))
        .interact()
        .context("flag predict")?;
    let dwell = Confirm::with_theme(theme)
        .with_prompt("  flag enable_dwell_time")
        .default(get_bool("enable_dwell_time", false))
        .interact()
        .context("flag enable_dwell_time")?;

    let mut map = match existing {
        Some(serde_yaml::Value::Mapping(m)) => m.clone(),
        _ => serde_yaml::Mapping::new(),
    };
    let set = |m: &mut serde_yaml::Mapping, k: &str, v: bool| {
        m.insert(serde_yaml::Value::String(k.to_string()), serde_yaml::Value::Bool(v));
    };
    set(&mut map, "add_indel_length", add_indel);
    set(&mut map, "predict", predict);
    set(&mut map, "enable_dwell_time", dwell);
    Ok(serde_yaml::Value::Mapping(map))
}

#[derive(Debug, Clone)]
struct RegistrySummary {
    entries: Vec<RegistryEntrySummary>,
}

#[derive(Debug, Clone)]
struct RegistryEntrySummary {
    model_id: String,
    ts_path: PathBuf,
    status: EntryStatus,
}

#[derive(Debug, Clone)]
enum EntryStatus {
    Ok { size_bytes: u64 },
    ShaMismatch { pinned: String },
    FileMissing,
    Error(String),
}

fn inspect_registry(yaml: &Path) -> RegistrySummary {
    let reg = match Registry::load(yaml) {
        Ok(r) => r,
        Err(e) => {
            return RegistrySummary {
                entries: vec![RegistryEntrySummary {
                    model_id: "<parse error>".to_string(),
                    ts_path: yaml.to_path_buf(),
                    status: EntryStatus::Error(format!("{e:#}")),
                }],
            };
        }
    };
    let yaml_dir = yaml.parent().unwrap_or_else(|| Path::new("."));
    let mut entries = Vec::with_capacity(reg.models.len());
    for entry in reg.models.iter() {
        let ts_path = yaml_dir.join(&entry.torchscript.path);
        let status = if !ts_path.exists() {
            EntryStatus::FileMissing
        } else {
            match reg.verify_torchscript(entry, yaml_dir) {
                Ok(_) => {
                    let size = std::fs::metadata(&ts_path).map(|m| m.len()).unwrap_or(0);
                    EntryStatus::Ok { size_bytes: size }
                }
                Err(e) => {
                    let msg = format!("{e:#}");
                    if msg.contains("sha256 mismatch") {
                        EntryStatus::ShaMismatch {
                            pinned: entry.torchscript.sha256.clone(),
                        }
                    } else {
                        EntryStatus::Error(msg)
                    }
                }
            }
        };
        entries.push(RegistryEntrySummary {
            model_id: entry.model_id.clone(),
            ts_path,
            status,
        });
    }
    RegistrySummary { entries }
}

fn print_registry_summary(s: &RegistrySummary) {
    let n_total = s.entries.len();
    let n_ok = s
        .entries
        .iter()
        .filter(|e| matches!(e.status, EntryStatus::Ok { .. }))
        .count();
    let n_missing = s
        .entries
        .iter()
        .filter(|e| matches!(e.status, EntryStatus::FileMissing))
        .count();
    let n_bad = n_total - n_ok - n_missing;
    println!(
        "  {} {} loadable, {} template (file missing), {} broken",
        style(format!("Registry: {n_total} entries")).bold(),
        style(n_ok).green(),
        style(n_missing).yellow(),
        style(n_bad).red(),
    );
    for e in &s.entries {
        let mark = match &e.status {
            EntryStatus::Ok { size_bytes } => style(format!(
                "  ✓ {} ({:.1} MB, sha OK)",
                e.model_id,
                *size_bytes as f64 / 1_048_576.0,
            ))
            .green(),
            EntryStatus::FileMissing => style(format!(
                "  · {} (template — file missing: {})",
                e.model_id,
                e.ts_path.display()
            ))
            .yellow(),
            EntryStatus::ShaMismatch { pinned } => style(format!(
                "  ✗ {} (sha mismatch; pinned: {})",
                e.model_id,
                &pinned[..pinned.len().min(16)]
            ))
            .red(),
            EntryStatus::Error(msg) => style(format!(
                "  ✗ {} (error: {})",
                e.model_id,
                &msg[..msg.len().min(80)]
            ))
            .red(),
        };
        println!("{mark}");
    }
}

fn print_entries_verbose(reg: &Registry) {
    for (i, entry) in reg.models.iter().enumerate() {
        println!();
        println!(
            "  {} {} ({})",
            style(format!("[{}]", i + 1)).dim(),
            style(&entry.model_id).bold(),
            entry.platform,
        );
        println!("      stage         {}", entry.stage);
        println!(
            "      input         shape={:?} dtype={}",
            entry.input.shape, entry.input.dtype
        );
        println!(
            "      output        shape={:?} dtype={}",
            entry.output.shape, entry.output.dtype
        );
        println!(
            "      torchscript   {} (sha {}…)",
            entry.torchscript.path,
            &entry.torchscript.sha256[..entry.torchscript.sha256.len().min(12)]
        );
    }
    println!();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_ghost_registry(dir: &Path) -> PathBuf {
        let yaml = dir.join("registry.yaml");
        std::fs::write(
            &yaml,
            "server_version: \"0.0.0\"\nmodels:\n  - model_id: ghost\n    platform: ont\n    stage: pileup\n    torchscript:\n      path: ghost.ts\n      sha256: 0000000000000000000000000000000000000000000000000000000000000000\n    input:\n      shape: [-1, 33, 18]\n      dtype: int32\n    output:\n      shape: [-1, 24]\n      dtype: float32\n    flags:\n      add_indel_length: false\n      predict: true\n",
        )
        .unwrap();
        yaml
    }

    fn tmp() -> PathBuf {
        let d = std::env::temp_dir().join(format!(
            "clair3-connect-radmin-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn inspect_marks_missing_file_as_template() {
        let dir = tmp();
        let yaml = write_ghost_registry(&dir);
        let s = inspect_registry(&yaml);
        assert_eq!(s.entries.len(), 1);
        assert!(matches!(s.entries[0].status, EntryStatus::FileMissing));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn add_save_reload_roundtrip() {
        let dir = tmp();
        let yaml = write_ghost_registry(&dir);
        let mut reg = Registry::load(&yaml).unwrap();
        let mut flags = serde_yaml::Mapping::new();
        flags.insert("predict".into(), true.into());
        reg.models.push(ModelEntry {
            model_id: "added_model".to_string(),
            description: String::new(),
            platform: "hifi".to_string(),
            stage: "full_alignment".to_string(),
            accuracy: Some("revio".to_string()),
            checkpoint_source: None,
            torchscript: TorchScript {
                path: "added.ts".to_string(),
                sha256: "1".repeat(64),
            },
            input: TensorSpec { shape: vec![-1, 55, 33, 8], dtype: "int8".to_string(), layout: None },
            output: TensorSpec { shape: vec![-1, 90], dtype: "float32".to_string(), layout: None },
            flags: serde_yaml::Value::Mapping(flags),
            notes: None,
        });
        save_registry(&reg, &yaml).unwrap();
        assert!(yaml.with_extension("yaml.bak").exists(), "backup not written");
        let reloaded = Registry::load(&yaml).unwrap();
        assert_eq!(reloaded.models.len(), 2);
        let added = reloaded.get("added_model").expect("added entry survived");
        assert_eq!(added.platform, "hifi");
        assert_eq!(added.input.shape, vec![-1, 55, 33, 8]);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn parse_shape_handles_spaces_and_negatives() {
        assert_eq!(parse_shape("-1, 89, 33, 8").unwrap(), vec![-1, 89, 33, 8]);
        assert_eq!(parse_shape("-1,24").unwrap(), vec![-1, 24]);
        assert!(parse_shape("").is_err());
        assert!(parse_shape("3,x").is_err());
    }
}
