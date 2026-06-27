#!/usr/bin/env python3
"""E4 APOE task definitions.

Each task has a natural-language `prompt`, a `match` rule, and an `expected_fn`
that returns the deterministic ground-truth answer. For APOE the truth is the
EXTERNAL 1000G phased-panel diplotype (a constant), so no Clair3-connect computation is
needed to freeze it — tasks.json is just these two tasks with their expected
diplotypes.
"""

def _apoe_task(sample, diplotype):
    """APOE ε-diplotyping task for a 1000G sample. Genotype-based — NO phasing:
    with the three standard haplotypes ε2/ε3/ε4, the two SNP genotypes map
    uniquely to a diplotype (a C at rs429358 can only sit on ε4; a T at rs7412
    only on ε2), so het+het is unambiguously ε2/ε4. Ground truth is the EXTERNAL
    1000G phased-panel diplotype (frozen), not Clair3-connect's own output."""
    return {
        "id": f"apoe_diplotype_{sample.lower()}",
        "prompt": (
            f"Clinical genomics. Sample {sample} is a human whole-genome sequencing "
            "dataset registered with the service. Determine its APOE diplotype — the "
            "ε2/ε3/ε4 combination that stratifies Alzheimer's-disease risk — defined "
            "by two SNPs on chr19 (GRCh38): rs429358 at chr19:44908684 (ref T, alt C) "
            "and rs7412 at chr19:44908822 (ref C, alt T). The three common haplotypes "
            "are ε4=(rs429358 C, rs7412 C), ε3=(T, C), ε2=(T, T). Steps: (1) determine "
            "which sequencing platform this sample is and select the matching Clair3 "
            "model; (2) call variants and read the genotype at BOTH SNP sites — a site "
            "with NO called variant is homozygous reference; (3) combine the two "
            "genotypes into the standard clinical diplotype. Finish with one line "
            "'ANSWER: <diplotype>', for example 'ANSWER: e4/e4'."
        ),
        "match": "diplotype",
        "expected_fn": (lambda d=diplotype: d),
    }

TASKS = [
    _apoe_task("HG00136", "e4/e4"),
    _apoe_task("HG00150", "e2/e4"),
]
