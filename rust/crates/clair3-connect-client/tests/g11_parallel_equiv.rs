use clair3_connect_client::chunk_planner::{plan_chunks, ChunkRegion};
use clair3_connect_client::parallel_driver::{run_parallel_pileup, ParallelPileupConfig};
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
async fn parallel_pileup_matches_single_threaded() {
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

    let mut single_buf: Vec<u8> = Vec::new();
    let single_inputs = PileupPipelineInputs {
        bam_path: &bam,
        fasta_path: &fa,
        region: expand_user_region("chr20", 100_000, 300_000),
        server_endpoint: endpoint.clone(),
        model_id: "pileup_ont_hac_v500".to_string(),
        sample_name: "HG003".to_string(),
        clair3_version: "2.0.1".to_string(),
        max_depth: 144,
        pileup_params: pileup_params.clone(),
        output_config,
        batch_size: 32,
        aead_key: None,
        tls: None,
    };
    run_pileup_pipeline(&single_inputs, &mut single_buf)
        .await
        .expect("single pileup");
    let single_text = String::from_utf8(single_buf).unwrap();
    let single_rows: Vec<&str> = single_text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect();
    eprintln!("[single] {} rows", single_rows.len());

    let mut chunks: Vec<ChunkRegion> = Vec::new();
    let chunk_size = 50_000_u64;
    let start = 100_000_u64;
    let end = 300_000_u64;
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

    assert_eq!(chunks.len(), 4);

    let _ = plan_chunks;

    let mut par_buf: Vec<u8> = Vec::new();
    let par_config = ParallelPileupConfig {
        bam_path: bam.clone(),
        fasta_path: fa.clone(),
        server_endpoint: endpoint.clone(),
        model_id: "pileup_ont_hac_v500".to_string(),
        sample_name: "HG003".to_string(),
        clair3_version: "2.0.1".to_string(),
        max_depth: 144,
        pileup_params,
        output_config,
        batch_size: 32,
        max_concurrent: 4,
        manifest_dir: None,
        verify_resume: false,
        aead_key: None,
        tls: None,
    };
    let stats = run_parallel_pileup(&chunks, &par_config, &mut par_buf)
        .await
        .expect("parallel pileup");
    let par_text = String::from_utf8(par_buf).unwrap();
    let par_rows: Vec<&str> = par_text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect();
    eprintln!(
        "[parallel] {} rows from {}/{} chunks ({} failed)",
        par_rows.len(),
        stats.chunks_completed,
        stats.chunks_total,
        stats.chunks_failed,
    );
    assert_eq!(stats.chunks_failed, 0);

    use std::collections::HashMap;
    let single_map: HashMap<&str, &str> = single_rows
        .iter()
        .map(|r| (r.split('\t').nth(1).unwrap_or(""), *r))
        .collect();
    let mut par_map: HashMap<&str, &str> = HashMap::new();
    let mut par_dups: Vec<&str> = Vec::new();
    for r in &par_rows {
        let pos = r.split('\t').nth(1).unwrap_or("");
        if par_map.insert(pos, *r).is_some() {
            par_dups.push(pos);
        }
    }
    eprintln!("[parallel] duplicate positions: {}", par_dups.len());
    assert!(
        par_dups.is_empty(),
        "parallel output has {} duplicate positions; first: {:?}",
        par_dups.len(),
        par_dups.first(),
    );

    let only_in_single: Vec<&str> = single_map
        .keys()
        .filter(|k| !par_map.contains_key(*k))
        .copied()
        .collect();
    let only_in_par: Vec<&str> = par_map
        .keys()
        .filter(|k| !single_map.contains_key(*k))
        .copied()
        .collect();
    eprintln!(
        "[diff] single-only: {} positions; parallel-only: {} positions",
        only_in_single.len(),
        only_in_par.len()
    );
    assert!(
        only_in_single.len() <= 5 && only_in_par.len() <= 5,
        "position sets differ too much"
    );

    let mut row_diffs = 0;
    let mut first_row_diff: Option<(String, String, String)> = None;
    for (pos, single_row) in &single_map {
        let Some(par_row) = par_map.get(pos) else {
            continue;
        };
        if single_row != par_row {
            row_diffs += 1;
            if first_row_diff.is_none() {
                first_row_diff =
                    Some((pos.to_string(), single_row.to_string(), par_row.to_string()));
            }
        }
    }
    if let Some((p, s, par)) = &first_row_diff {
        eprintln!("[first row diff @ pos {}]\n  single:   {}\n  parallel: {}", p, s, par);
    }
    let overlap = single_map.len().min(par_map.len());
    let row_match = overlap - row_diffs;
    let ratio = row_match as f64 / single_map.len() as f64;
    eprintln!(
        "[row-parity] {}/{} byte-equal ({:.3}%); {} row diffs",
        row_match,
        single_map.len(),
        100.0 * ratio,
        row_diffs,
    );
    assert!(
        ratio >= 0.95,
        "row-level parity below 95%: {}/{}",
        row_match,
        single_map.len()
    );

    server_handle.abort();
}
