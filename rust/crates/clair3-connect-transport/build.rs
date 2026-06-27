fn main() -> Result<(), Box<dyn std::error::Error>> {
    let proto_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .ok_or("missing parent")?
        .parent()
        .ok_or("missing grandparent")?
        .join("proto")
        .join("clair3-connect.proto");
    let proto_dir = proto_path
        .parent()
        .ok_or("missing proto dir")?
        .to_path_buf();

    tonic_build::configure()
        .build_server(true)
        .build_client(true)

        .bytes([
            ".clair3_connect.v1.PredictRequest.payload",
            ".clair3_connect.v1.PredictResponse.y_payload",
            ".clair3_connect.v1.PredictRequest.aead_nonce",
            ".clair3_connect.v1.PredictResponse.y_aead_nonce",
        ])
        .compile_protos(&[&proto_path], &[&proto_dir])?;
    println!("cargo:rerun-if-changed={}", proto_path.display());
    Ok(())
}
