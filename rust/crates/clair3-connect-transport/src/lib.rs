#![deny(rust_2018_idioms)]
#![warn(missing_docs)]

#[allow(missing_docs)]
pub mod proto {
    tonic::include_proto!("clair3_connect.v1");
}

pub mod aead;
pub mod checksum;
pub mod codec;
pub mod predict_client;
pub mod tls;
