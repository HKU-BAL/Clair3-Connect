use anyhow::{anyhow, Context, Result};
use std::path::Path;
use tch::{CModule, Kind, Tensor};

use crate::model_registry::{ModelEntry, TensorSpec};

pub struct LoadedModel {

    pub module: CModule,

    pub entry: ModelEntry,

    pub device: tch::Device,
}

impl LoadedModel {

    pub fn load(ts_path: &Path, entry: ModelEntry, device: tch::Device) -> Result<Self> {
        let mut module = CModule::load_on_device(ts_path, device)
            .with_context(|| format!("load CModule {}", ts_path.display()))?;
        module.set_eval();
        Ok(Self {
            module,
            entry,
            device,
        })
    }

    pub fn validate_input(&self, x: &Tensor) -> Result<()> {
        validate_against_spec(&self.entry.model_id, &self.entry.input, &x.size(), x.kind())
    }

    pub fn forward(&self, x: &Tensor) -> Result<Tensor> {
        self.validate_input(x)?;
        let y = tch::no_grad(|| self.module.forward_ts(&[x]))
            .with_context(|| format!("forward {}", self.entry.model_id))?;
        Ok(y)
    }
}

pub fn parse_dtype(s: &str) -> Result<Kind> {
    Ok(match s {
        "int8" => Kind::Int8,
        "uint8" => Kind::Uint8,
        "int16" => Kind::Int16,
        "int32" => Kind::Int,
        "int64" => Kind::Int64,
        "float16" => Kind::Half,
        "float32" => Kind::Float,
        "float64" => Kind::Double,
        other => return Err(anyhow!("unknown dtype {other:?}")),
    })
}

pub fn validate_against_spec(
    model_id: &str,
    spec: &TensorSpec,
    actual_shape: &[i64],
    actual_dtype: Kind,
) -> Result<()> {
    if actual_shape.len() != spec.shape.len() {
        return Err(anyhow!(
            "model {}: input rank mismatch (expected {} dims, got {})",
            model_id,
            spec.shape.len(),
            actual_shape.len()
        ));
    }
    for (i, (&exp, &act)) in spec.shape.iter().zip(actual_shape.iter()).enumerate() {
        if i == 0 && exp == -1 {
            if act <= 0 {
                return Err(anyhow!(
                    "model {}: batch dim must be positive (got {})",
                    model_id,
                    act
                ));
            }
            continue;
        }
        if exp == -1 {
            return Err(anyhow!(
                "model {}: registry declares -1 at non-batch dim {} ({:?}); \
                 only dim 0 (batch) may be dynamic",
                model_id,
                i,
                spec.shape
            ));
        }
        if exp != act {
            return Err(anyhow!(
                "model {}: shape mismatch at dim {} (expected {}, got {})",
                model_id,
                i,
                exp,
                act
            ));
        }
    }
    let expected_kind = parse_dtype(&spec.dtype)?;
    if actual_dtype != expected_kind {
        return Err(anyhow!(
            "model {}: dtype mismatch (expected {} / {:?}, got {:?})",
            model_id,
            spec.dtype,
            expected_kind,
            actual_dtype
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model_registry::TensorSpec;

    fn spec(shape: Vec<i64>, dtype: &str) -> TensorSpec {
        TensorSpec {
            shape,
            dtype: dtype.to_string(),
            layout: None,
        }
    }

    #[test]
    fn accepts_dynamic_batch() {
        let s = spec(vec![-1, 33, 18], "int32");
        validate_against_spec("m", &s, &[988, 33, 18], Kind::Int).unwrap();
        validate_against_spec("m", &s, &[1, 33, 18], Kind::Int).unwrap();
        validate_against_spec("m", &s, &[10_000, 33, 18], Kind::Int).unwrap();
    }

    #[test]
    fn rejects_non_batch_negative_one() {
        let s = spec(vec![-1, -1, 18], "int32");
        let err = validate_against_spec("m", &s, &[10, 33, 18], Kind::Int).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("non-batch dim 1"),
            "expected non-batch-dim error, got: {msg}"
        );
    }

    #[test]
    fn rejects_rank_mismatch() {
        let s = spec(vec![-1, 33, 18], "int32");
        let err = validate_against_spec("m", &s, &[10, 33], Kind::Int).unwrap_err();
        assert!(err.to_string().contains("rank mismatch"));
    }

    #[test]
    fn rejects_shape_mismatch() {
        let s = spec(vec![-1, 33, 18], "int32");
        let err = validate_against_spec("m", &s, &[10, 32, 18], Kind::Int).unwrap_err();
        assert!(err.to_string().contains("shape mismatch at dim 1"));
    }

    #[test]
    fn rejects_dtype_mismatch() {
        let s = spec(vec![-1, 33, 18], "int32");
        let err = validate_against_spec("m", &s, &[10, 33, 18], Kind::Float).unwrap_err();
        assert!(err.to_string().contains("dtype mismatch"));
    }

    #[test]
    fn rejects_zero_or_negative_batch() {
        let s = spec(vec![-1, 33, 18], "int32");
        assert!(validate_against_spec("m", &s, &[0, 33, 18], Kind::Int).is_err());

        assert!(validate_against_spec("m", &s, &[-1, 33, 18], Kind::Int).is_err());
    }

    #[test]
    fn fixed_batch_must_match() {
        let s = spec(vec![4, 33, 18], "int32");
        validate_against_spec("m", &s, &[4, 33, 18], Kind::Int).unwrap();
        assert!(validate_against_spec("m", &s, &[8, 33, 18], Kind::Int).is_err());
    }

    #[test]
    fn fa_shape_int8() {
        let s = spec(vec![-1, 89, 33, 8], "int8");
        validate_against_spec("fa", &s, &[16, 89, 33, 8], Kind::Int8).unwrap();
        assert!(validate_against_spec("fa", &s, &[16, 89, 33, 9], Kind::Int8).is_err());
    }
}
