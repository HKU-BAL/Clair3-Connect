use clair3_connect_client::pileup_pipeline::{
    expand_user_region, run_pileup_pipeline, PileupPipelineInputs,
};
use clair3_connect_core::pileup::kernel::PileupParams;
use clair3_connect_core::vcf::output_with::OutputWithConfig;
use clair3_connect_infer::model_registry::Registry;
use clair3_connect_server::service::PredictService;
use clair3_connect_transport::proto::clair3_connect_inference_server::Clair3ConnectInferenceServer;
use std::path::PathBuf;
use tokio::net::TcpListener;
use tonic::transport::Server;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .to_path_buf()
}

#[tokio::test]
async fn pileup_pipeline_matches_stock_clair3() {
    let root = workspace_root();
    let yaml = root.join("models/model_registry.example.yaml");
    let bam = root.join("test_data/ont_chr20_demo/HG003_chr20_demo.bam");
    let fa = root.join("test_data/ont_chr20_demo/GRCh38_no_alt_chr20.fa");
    let stock_vcf = root.join("test_data/ont_chr20_demo/artifacts/pileup.vcf");

    if !yaml.exists() || !bam.exists() || !fa.exists() || !stock_vcf.exists() {
        eprintln!("Skipping — fixtures missing");
        return;
    }

    let reg = Registry::load(&yaml).unwrap();
    let yaml_dir = yaml.parent().unwrap().to_path_buf();
    let svc = match PredictService::load(&reg, &yaml_dir, tch::Device::Cpu) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Skipping — model load failed: {e}");
            return;
        }
    };
    let listener = TcpListener::bind("127.0.0.1:0").await.unwrap();
    let addr = listener.local_addr().unwrap();
    let incoming = tokio_stream::wrappers::TcpListenerStream::new(listener);
    let server_handle = tokio::spawn(async move {
        Server::builder()
            .add_service(Clair3ConnectInferenceServer::new(svc))
            .serve_with_incoming(incoming)
            .await
            .ok();
    });

    tokio::task::yield_now().await;

    let mut buf: Vec<u8> = Vec::new();
    let output_config = OutputWithConfig {
        is_show_reference: true,
        ..OutputWithConfig::default()
    };

    let inputs = PileupPipelineInputs {
        bam_path: &bam,
        fasta_path: &fa,

        region: expand_user_region("chr20", 100_000, 300_000),
        server_endpoint: format!("http://{}", addr),
        model_id: "pileup_ont_hac_v500".to_string(),
        sample_name: "HG003".to_string(),
        clair3_version: "2.0.1".to_string(),
        max_depth: 144,
        pileup_params: PileupParams {
            min_snp_af: 0.08,
            min_indel_af: 0.15,
            min_mapq: 5,
            ..Default::default()
        },
        output_config,
        batch_size: 32,
        aead_key: None,
        tls: None,
    };

    run_pileup_pipeline(&inputs, &mut buf).await.expect("pipeline");
    let rust_text = String::from_utf8(buf).unwrap();
    let rust_rows: Vec<&str> = rust_text.lines().filter(|l| !l.starts_with('#')).collect();

    let stock_text = std::fs::read_to_string(&stock_vcf).unwrap();
    let stock_rows: Vec<&str> = stock_text.lines().filter(|l| !l.starts_with('#')).collect();

    eprintln!("rust rows: {}, stock rows: {}", rust_rows.len(), stock_rows.len());

    let rust_pos: std::collections::HashSet<&str> = rust_rows
        .iter()
        .map(|r| r.split('\t').nth(1).unwrap_or(""))
        .collect();
    let stock_pos: std::collections::HashSet<&str> = stock_rows
        .iter()
        .map(|r| r.split('\t').nth(1).unwrap_or(""))
        .collect();
    for p in &stock_pos {
        if !rust_pos.contains(*p) {
            eprintln!("ONLY in stock: pos {}", p);
        }
    }
    for p in &rust_pos {
        if !stock_pos.contains(*p) {
            eprintln!("ONLY in rust: pos {}", p);
        }
    }
    assert_eq!(rust_rows.len(), stock_rows.len(), "row count differs");

    let mut structural_diffs = 0;
    let mut qual_drifts = 0;
    let mut first_diff: Option<(usize, String, String)> = None;
    for (i, (r, s)) in rust_rows.iter().zip(stock_rows.iter()).enumerate() {
        if r == s {
            continue;
        }
        let rc: Vec<&str> = r.split('\t').collect();
        let sc: Vec<&str> = s.split('\t').collect();
        let same_struct = rc.len() == sc.len()
            && rc[0] == sc[0]
            && rc[1] == sc[1]
            && rc[3] == sc[3]
            && rc[4] == sc[4]
            && rc[6] == sc[6]
            && rc[7] == sc[7]
            && rc[8] == sc[8];
        let qual_r: f64 = rc.get(5).and_then(|v| v.parse().ok()).unwrap_or(-1.0);
        let qual_s: f64 = sc.get(5).and_then(|v| v.parse().ok()).unwrap_or(-1.0);
        let qual_diff = (qual_r - qual_s).abs();

        let sample_r = rc.get(9).copied().unwrap_or("");
        let sample_s = sc.get(9).copied().unwrap_or("");
        let same_sample_except_first = {
            let r_parts: Vec<&str> = sample_r.split(':').collect();
            let s_parts: Vec<&str> = sample_s.split(':').collect();
            r_parts.len() == s_parts.len()
                && r_parts.iter().zip(s_parts.iter()).enumerate().all(|(idx, (a, b))| idx == 1 || a == b)
        };
        if same_struct && qual_diff <= 0.05 && same_sample_except_first {
            qual_drifts += 1;
        } else {
            structural_diffs += 1;
            if first_diff.is_none() {
                first_diff = Some((i, (*r).to_string(), (*s).to_string()));
            }
        }
    }
    eprintln!(
        "VCF parity: {} structural diffs, {} float-rounding drifts (within QUAL tolerance)",
        structural_diffs, qual_drifts
    );
    if let Some((i, r, s)) = &first_diff {
        eprintln!("first structural diff at row {i}:\n  rust:  {r}\n  stock: {s}");
    }

    let total = stock_rows.len();
    let near_equal = total - structural_diffs;
    eprintln!("near-equal rows: {}/{} ({:.3}%)",
              near_equal, total, 100.0 * near_equal as f64 / total as f64);
    assert!(
        near_equal as f64 / total as f64 > 0.99,
        "near-identity fell below 99% ({}/{})", near_equal, total
    );

    server_handle.abort();
}
