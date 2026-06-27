#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Platform {

    Ont,

    Hifi,

    Ilmn,
}

impl Platform {

    pub fn from_str_or_ont(s: &str) -> Self {
        match s {
            "ont" => Platform::Ont,
            "hifi" => Platform::Hifi,
            "ilmn" | "illumina" => Platform::Ilmn,
            other => {
                tracing::warn!(
                    platform = %other,
                    "unknown platform; defaulting to ONT"
                );
                Platform::Ont
            }
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Platform::Ont => "ont",
            Platform::Hifi => "hifi",
            Platform::Ilmn => "ilmn",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct PlatformParams {

    pub max_depth: usize,

    pub matrix_depth: usize,

    pub min_af: f64,

    pub min_mq: u8,

    pub longphase_platform_flag: &'static str,
}

impl PlatformParams {

    pub const fn for_platform(platform: Platform) -> Self {
        match platform {
            Platform::Ont => Self {
                max_depth: 144,
                matrix_depth: 89,
                min_af: 0.15,
                min_mq: 5,
                longphase_platform_flag: "--ont",
            },
            Platform::Hifi => Self {
                max_depth: 144,
                matrix_depth: 55,
                min_af: 0.08,
                min_mq: 5,
                longphase_platform_flag: "--pb",
            },
            Platform::Ilmn => Self {
                max_depth: 144,
                matrix_depth: 55,
                min_af: 0.08,
                min_mq: 5,

                longphase_platform_flag: "--ont",
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ont_params() {
        let p = PlatformParams::for_platform(Platform::Ont);
        assert_eq!(p.max_depth, 144);
        assert_eq!(p.matrix_depth, 89);
        assert_eq!(p.min_af, 0.15);
        assert_eq!(p.longphase_platform_flag, "--ont");
    }

    #[test]
    fn hifi_params() {
        let p = PlatformParams::for_platform(Platform::Hifi);
        assert_eq!(p.max_depth, 144);
        assert_eq!(p.matrix_depth, 55);
        assert_eq!(p.min_af, 0.08);
        assert_eq!(p.longphase_platform_flag, "--pb");
    }

    #[test]
    fn ilmn_params() {
        let p = PlatformParams::for_platform(Platform::Ilmn);
        assert_eq!(p.max_depth, 144);
        assert_eq!(p.matrix_depth, 55);
        assert_eq!(p.min_af, 0.08);
    }

    #[test]
    fn parse_known_names() {
        assert_eq!(Platform::from_str_or_ont("ont"), Platform::Ont);
        assert_eq!(Platform::from_str_or_ont("hifi"), Platform::Hifi);
        assert_eq!(Platform::from_str_or_ont("ilmn"), Platform::Ilmn);
        assert_eq!(Platform::from_str_or_ont("illumina"), Platform::Ilmn);
    }

    #[test]
    fn parse_unknown_defaults_to_ont() {
        assert_eq!(Platform::from_str_or_ont("rubbish"), Platform::Ont);
    }

    #[test]
    fn as_str_roundtrip() {
        for p in [Platform::Ont, Platform::Hifi, Platform::Ilmn] {
            assert_eq!(Platform::from_str_or_ont(p.as_str()), p);
        }
    }
}
