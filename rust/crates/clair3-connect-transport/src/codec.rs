use crate::proto::TensorDType;
use std::mem::size_of;

pub fn expected_bytes(shape: &[i64], dtype: TensorDType) -> Option<usize> {
    let elt: usize = match dtype {
        TensorDType::TensorDtypeInt8 => size_of::<i8>(),
        TensorDType::TensorDtypeInt16 => size_of::<i16>(),
        TensorDType::TensorDtypeInt32 => size_of::<i32>(),
        TensorDType::TensorDtypeFloat32 => size_of::<f32>(),
        TensorDType::TensorDtypeUnspecified => return None,
    };
    let mut n: usize = 1;
    for d in shape {
        if *d < 0 {
            return None;
        }
        n = n.checked_mul(*d as usize)?;
    }
    n.checked_mul(elt)
}

pub fn encode_i8(values: &[i8]) -> Vec<u8> {
    values.iter().map(|v| *v as u8).collect()
}

pub fn encode_i32(values: &[i32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(values.len() * 4);
    for v in values {
        out.extend_from_slice(&v.to_le_bytes());
    }
    out
}

pub fn encode_f32(values: &[f32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(values.len() * 4);
    for v in values {
        out.extend_from_slice(&v.to_le_bytes());
    }
    out
}

pub fn decode_f32(bytes: &[u8]) -> Vec<f32> {
    let mut out = Vec::with_capacity(bytes.len() / 4);
    for chunk in bytes.chunks_exact(4) {
        out.push(f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]));
    }
    out
}

pub fn decode_i32(bytes: &[u8]) -> Vec<i32> {
    let mut out = Vec::with_capacity(bytes.len() / 4);
    for chunk in bytes.chunks_exact(4) {
        out.push(i32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]));
    }
    out
}

pub fn decode_i8(bytes: &[u8]) -> Vec<i8> {
    bytes.iter().map(|b| *b as i8).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expected_bytes_simple() {
        assert_eq!(
            expected_bytes(&[1, 33, 18], TensorDType::TensorDtypeInt32),
            Some(33 * 18 * 4)
        );
        assert_eq!(
            expected_bytes(&[4, 89, 33, 8], TensorDType::TensorDtypeInt8),
            Some(4 * 89 * 33 * 8)
        );
    }

    #[test]
    fn round_trip_f32() {
        let v: Vec<f32> = vec![0.0, -1.5, 9.5_f32, f32::INFINITY, f32::NAN, -0.0];
        let bytes = encode_f32(&v);
        let back = decode_f32(&bytes);
        assert_eq!(back.len(), v.len());
        for (a, b) in v.iter().zip(back.iter()) {
            if a.is_nan() {
                assert!(b.is_nan());
            } else {
                assert_eq!(a.to_bits(), b.to_bits());
            }
        }
    }

    #[test]
    fn round_trip_i32_i8() {
        let i32s: Vec<i32> = vec![-5, 0, 1000, i32::MIN, i32::MAX];
        assert_eq!(decode_i32(&encode_i32(&i32s)), i32s);
        let i8s: Vec<i8> = vec![-128, -1, 0, 1, 127];
        assert_eq!(decode_i8(&encode_i8(&i8s)), i8s);
    }
}

pub const PAYLOAD_ENCODING_ZSTD: &str = "zstd";

pub fn request_zstd_enabled() -> bool {
    !matches!(
        std::env::var("CLAIR3_CONNECT_COMPRESS").ok().as_deref(),
        Some("none") | Some("off") | Some("0") | Some("gzip")
    )
}

pub fn compress_payload(payload: &[u8]) -> std::io::Result<Vec<u8>> {
    zstd::bulk::compress(payload, 1)
}

pub fn decompress_payload(payload: &[u8], max_len: usize) -> std::io::Result<Vec<u8>> {
    zstd::bulk::decompress(payload, max_len)
}
