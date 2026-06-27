use clair3_connect_client::chunk_planner::ChunkRegion;
use clair3_connect_client::parallel_driver::{run_parallel_pileup, ParallelPileupConfig};
use clair3_connect_core::manifest::read_manifest;
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

fn build_chunks() -> Vec<ChunkRegion> {
    let mut chunks: Vec<ChunkRegion> = Vec::new();
    let start = 100_000_u64;
    let end = 300_000_u64;
    let chunk_size = 50_000_u64;
    let chunk_count = ((end - start) / chunk_size) as u32;
    for i in 0..chunk_count {
        let cs = start + (i as u64) * chunk_size;
        let ce = (cs + chunk_size - 1).min(end);
        chunks.push(ChunkRegion {
            contig: "chr20".to_string(),
            start_1based: cs,
            end_1based: ce,
            chunk_idx: i,
            chunk_count,
        });
    }
    chunks
}

#[tokio::test]
async fn resume_yields_identical_vcf() {
    let root = workspace_root();
    let yaml = root.join("models/model_registry.example.yaml");
    let bam = root.join("test_data/ont_chr20_demo/HG003_chr20_demo.bam");
    let fa = root.join("test_data/ont_chr20_demo/GRCh38_no_alt_chr20.fa");
    if !yaml.exists() || !bam.exists() || !fa.exists() {
        eprintln!("Skipping — fixtures missing");
        return;
    }

    let reg = Registry::load(&yaml).unwrap();
    let yaml_dir = yaml.parent().unwrap().to_path_buf();
    let svc = match PredictService::load(&reg, &yaml_dir, tch::Device::Cpu) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Skipping — model load: {e}");
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
    let endpoint = format!("http://{}", addr);

    let manifest_tmp = tempfile::Builder::new()
        .prefix("clair3-connect-g12-resume-")
        .tempdir()
        .expect("create temp dir");
    let manifest_dir = manifest_tmp.path().to_path_buf();
    let chunks = build_chunks();
    let pileup_params = PileupParams {
        min_snp_af: 0.08,
        min_indel_af: 0.15,
        min_mapq: 5,
        ..Default::default()
    };
    let output_config = OutputWithConfig {
        is_show_reference: true,
        ..OutputWithConfig::default()
    };

    let mk_config = |dir: Option<PathBuf>| ParallelPileupConfig {
        bam_path: bam.clone(),
        fasta_path: fa.clone(),
        server_endpoint: endpoint.clone(),
        model_id: "pileup_ont_hac_v500".to_string(),
        sample_name: "HG003".to_string(),
        clair3_version: "2.0.1".to_string(),
        max_depth: 144,
        pileup_params: pileup_params.clone(),
        output_config,
        batch_size: 32,
        max_concurrent: 4,
        manifest_dir: dir,
        verify_resume: true,
        aead_key: None,
        tls: None,
    };

    let mut first_buf: Vec<u8> = Vec::new();
    let first_stats = run_parallel_pileup(
        &chunks,
        &mk_config(Some(manifest_dir.clone())),
        &mut first_buf,
    )
    .await
    .expect("first run");
    eprintln!(
        "[first] chunks={} completed={} resumed={} rows={}",
        first_stats.chunks_total,
        first_stats.chunks_completed,
        first_stats.chunks_resumed,
        first_stats.rows_emitted,
    );
    assert_eq!(first_stats.chunks_resumed, 0);
    assert_eq!(first_stats.chunks_failed, 0);

    let first_text = String::from_utf8(first_buf).unwrap();
    let first_rows: Vec<&str> = first_text.lines().collect();

    let entries = read_manifest(&manifest_dir).unwrap();
    assert_eq!(entries.len(), 4);
    let to_drop = &entries[2..];
    for e in to_drop {
        std::fs::remove_file(manifest_dir.join(&e.vcf_path)).ok();
    }

    let manifest_path = manifest_dir.join("manifest.jsonl");
    let text = std::fs::read_to_string(&manifest_path).unwrap();
    let kept_lines: Vec<&str> = text.lines().take(2).collect();
    let truncated = kept_lines.join("\n") + "\n";
    std::fs::write(&manifest_path, truncated).unwrap();
    eprintln!("[simulate kill] manifest truncated to 2/4 chunks");

    let mut second_buf: Vec<u8> = Vec::new();
    let second_stats = run_parallel_pileup(
        &chunks,
        &mk_config(Some(manifest_dir.clone())),
        &mut second_buf,
    )
    .await
    .expect("second run");
    eprintln!(
        "[second] chunks={} completed={} resumed={} rows={}",
        second_stats.chunks_total,
        second_stats.chunks_completed,
        second_stats.chunks_resumed,
        second_stats.rows_emitted,
    );
    assert_eq!(second_stats.chunks_resumed, 2, "should resume 2 chunks");
    assert_eq!(second_stats.chunks_completed, 4, "must complete all 4");

    let second_text = String::from_utf8(second_buf).unwrap();
    let second_rows: Vec<&str> = second_text.lines().collect();

    eprintln!(
        "[parity] first rows={} second rows={}",
        first_rows.len(),
        second_rows.len()
    );
    assert_eq!(
        first_rows.len(),
        second_rows.len(),
        "row count differs between original and resumed runs"
    );
    let mut diffs = 0;
    for (i, (a, b)) in first_rows.iter().zip(second_rows.iter()).enumerate() {
        if a != b {
            diffs += 1;
            if diffs <= 3 {
                eprintln!("row {i}:\n  first:  {a}\n  second: {b}");
            }
        }
    }
    assert_eq!(diffs, 0, "G12 failed: {diffs} rows differ");

    eprintln!("✅ G12 PASS: resumed VCF byte-identical to first run");

    drop(manifest_tmp);
    server_handle.abort();
}
