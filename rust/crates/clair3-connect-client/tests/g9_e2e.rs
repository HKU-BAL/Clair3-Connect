use clair3_connect_client::fa_pipeline::{run_fa_pipeline, FaPipelineInputs};
use clair3_connect_client::merge_vcf::{merge_vcf, MergeConfig};
use clair3_connect_client::pileup_pipeline::{
    expand_user_region, run_pileup_pipeline, PileupPipelineInputs,
};
use clair3_connect_client::select_candidates::{select_candidates, CandidateFilterMode};
use clair3_connect_client::select_qual::select_qual;
use clair3_connect_core::hts::region::Region;
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
async fn g9_full_mvp_pipeline() {
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

    let mut pileup_buf: Vec<u8> = Vec::new();
    let pileup_inputs = PileupPipelineInputs {
        bam_path: &bam,
        fasta_path: &fa,
        region: expand_user_region("chr20", 100_000, 300_000),
        server_endpoint: endpoint.clone(),
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
        output_config: OutputWithConfig {
            is_show_reference: true,
            ..OutputWithConfig::default()
        },
        batch_size: 32,
        aead_key: None,
        tls: None,
    };
    run_pileup_pipeline(&pileup_inputs, &mut pileup_buf).await.expect("pileup");
    let pileup_text = String::from_utf8(pileup_buf).unwrap();
    let pileup_rows: Vec<&str> = pileup_text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect();
    eprintln!("[1] pileup: {} rows", pileup_rows.len());
    assert!(pileup_rows.len() > 900, "pileup row count fell off");

    let qual_sel =
        select_qual(pileup_text.as_bytes(), 0.3, 0.3).expect("select_qual");
    eprintln!(
        "[2] select_qual: var_cutoff={:.2}, ref_cutoff={:.2}",
        qual_sel.var_qual_cut_off, qual_sel.ref_qual_cut_off
    );

    let cand_sel = select_candidates(
        pileup_text.as_bytes(),
        Some("chr20"),
        CandidateFilterMode::QualCutoffs {
            var_cutoff: qual_sel.var_qual_cut_off,
            ref_cutoff: qual_sel.ref_qual_cut_off,
        },
    )
    .expect("select_candidates");
    eprintln!(
        "[3] select_candidates: {} positions ({} variants, {} refs)",
        cand_sel.positions.len(),
        cand_sel.variants_kept,
        cand_sel.refs_kept,
    );
    assert!(!cand_sel.positions.is_empty(), "no FA candidates");

    let n_test = cand_sel.positions.len().min(64);
    let mut candidates_0b: Vec<u64> = cand_sel
        .positions
        .iter()
        .take(n_test)
        .map(|p| p - 1)
        .collect();
    candidates_0b.sort();
    candidates_0b.dedup();
    let mut fa_buf: Vec<u8> = Vec::new();
    let fa_inputs = FaPipelineInputs {
        bam_path: &bam,
        fasta_path: &fa,
        region: Region {
            contig: "chr20".to_string(),
            start: 99_999,
            end: 300_000,
        },
        candidates: candidates_0b.clone(),
        phased_variants: None,
        server_endpoint: endpoint.clone(),
        model_id: "full_alignment_ont_hac_v500".to_string(),
        sample_name: "HG003".to_string(),
        clair3_version: "2.0.1".to_string(),
        matrix_depth: 89,
        max_indel_length: 50,
        min_mq: 5,
        output_config: OutputWithConfig {
            is_show_reference: false,
            pileup: false,
            add_indel_length: true,
            ..OutputWithConfig::default()
        },
        batch_size: 32,
        aead_key: None,
        tls: None,
        low_mem: false,
        fa_window_bp: 0,
        streaming: false,
    };
    run_fa_pipeline(&fa_inputs, &mut fa_buf).await.expect("fa pipeline");
    let fa_text = String::from_utf8(fa_buf).unwrap();
    let fa_rows: Vec<&str> = fa_text.lines().filter(|l| !l.starts_with('#')).collect();
    eprintln!("[4] fa pipeline: {} rows from {} candidates", fa_rows.len(), candidates_0b.len());

    assert!(!fa_rows.is_empty(), "FA pipeline produced no rows");

    let mut merge_buf: Vec<u8> = Vec::new();
    let stats = merge_vcf(
        pileup_text.as_bytes(),
        fa_text.as_bytes(),
        &mut merge_buf,
        Some("chr20"),
        MergeConfig {
            quality_score_for_pass: Some(2.0),
            print_ref: false,
        },
    )
    .expect("merge");
    eprintln!(
        "[5] merge: {} fa rows + {} pileup rows = {} total (pileup overridden: {})",
        stats.fa_rows,
        stats.pileup_rows,
        stats.fa_rows + stats.pileup_rows,
        stats.pileup_overridden,
    );
    let merge_text = String::from_utf8(merge_buf).unwrap();
    let merge_rows: Vec<&str> = merge_text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect();
    eprintln!("[5] merge output: {} rows", merge_rows.len());

    assert_eq!(
        merge_rows.len(),
        stats.fa_rows + stats.pileup_rows,
        "merge row count off"
    );

    assert!(merge_rows.len() > 1, "merge_output is empty");

    server_handle.abort();
}
