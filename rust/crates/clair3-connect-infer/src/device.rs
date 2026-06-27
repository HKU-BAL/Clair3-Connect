use anyhow::{anyhow, Context, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Device {

    Cpu,

    Cuda(usize),
}

impl Device {

    pub fn parse(s: &str) -> Result<Self> {
        match s {
            "cpu" => Ok(Self::Cpu),
            "cuda" => Ok(Self::Cuda(0)),
            other => {
                if let Some(rest) = other.strip_prefix("cuda:") {
                    let n: usize = rest
                        .parse()
                        .with_context(|| format!("parse cuda index from {other:?}"))?;
                    Ok(Self::Cuda(n))
                } else {
                    Err(anyhow!(
                        "unknown device {other:?}; expected cpu | cuda | cuda:N"
                    ))
                }
            }
        }
    }

    pub fn as_tch(self) -> tch::Device {
        match self {
            Self::Cpu => tch::Device::Cpu,
            Self::Cuda(n) => tch::Device::Cuda(n),
        }
    }

    pub fn is_cuda(self) -> bool {
        matches!(self, Self::Cuda(_))
    }

    pub fn parse_list(s: &str) -> Result<Vec<Self>> {
        let mut out = Vec::new();
        for piece in s.split(',') {
            let t = piece.trim();
            if t.is_empty() {
                continue;
            }
            out.push(Self::parse(t)?);
        }
        if out.is_empty() {
            return Err(anyhow!("device list {s:?} parsed to zero devices"));
        }
        Ok(out)
    }
}

impl std::fmt::Display for Device {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Cpu => write!(f, "cpu"),
            Self::Cuda(n) => write!(f, "cuda:{n}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_cpu() {
        assert_eq!(Device::parse("cpu").unwrap(), Device::Cpu);
    }
    #[test]
    fn parse_cuda_default() {
        assert_eq!(Device::parse("cuda").unwrap(), Device::Cuda(0));
    }
    #[test]
    fn parse_cuda_n() {
        assert_eq!(Device::parse("cuda:3").unwrap(), Device::Cuda(3));
    }
    #[test]
    fn parse_invalid() {
        assert!(Device::parse("gpu").is_err());
        assert!(Device::parse("cuda:x").is_err());
    }
    #[test]
    fn display_roundtrip() {
        for s in ["cpu", "cuda:0", "cuda:7"] {
            assert_eq!(s, Device::parse(s).unwrap().to_string());
        }
    }

    #[test]
    fn parse_list_single() {
        assert_eq!(Device::parse_list("cuda:0").unwrap(), vec![Device::Cuda(0)]);
    }

    #[test]
    fn parse_list_multiple() {
        assert_eq!(
            Device::parse_list("cuda:0,cuda:1").unwrap(),
            vec![Device::Cuda(0), Device::Cuda(1)]
        );
    }

    #[test]
    fn parse_list_skips_empty() {
        assert_eq!(
            Device::parse_list("cuda:0,, cuda:1 ,").unwrap(),
            vec![Device::Cuda(0), Device::Cuda(1)]
        );
    }

    #[test]
    fn parse_list_empty_is_error() {
        assert!(Device::parse_list("").is_err());
        assert!(Device::parse_list(",,").is_err());
    }

    #[test]
    fn parse_list_bad_entry_is_error() {
        assert!(Device::parse_list("cuda:0,gpu").is_err());
    }
}
