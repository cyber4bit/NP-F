# ============================================================
# Snakemake Workflow: Network Pharmacology Pipeline
# ============================================================
# Example: 血筒 (Xuetong) → Rheumatoid Arthritis (RA)
#
# Usage:
#   python scripts/python/00_preflight.py --config config.yaml  # check first
#   snakemake --cores 8              # full pipeline
#   snakemake step01_compounds       # single step
#   snakemake --dry-run              # preview
# ============================================================

configfile: "config.yaml"

# ── Pipeline targets ──────────────────────────────────────────
DRUG    = config["drug"]["name"].replace(" ", "_")
DISEASE = config["disease"]["name"].replace(" ", "_")

rule all:
    input:
        "data/processed/01_compounds_filtered.csv",
        "data/processed/02_targets_merged.csv",
        "data/processed/03_disease_targets.csv",
        "data/processed/04_intersection_genes.csv",
        "data/processed/04_ppi_edges.csv",
        "results/figures/04_venn.pdf",
        "data/processed/05_hub_genes.csv",
        "results/figures/06_enrichment_kegg.pdf",
        "data/processed/07_mr_results.csv",
        "results/figures/07_mr_forest.pdf",
        "data/processed/08_docking_scores.csv",
        "results/figures/10_summary_figure.pdf"


# ── Step 00: Preflight check ────────────────────────────────
rule step00_preflight:
    output: touch("logs/00_preflight.done")
    log:    "logs/00_preflight.log"
    shell:
        "python scripts/python/00_preflight.py --config config.yaml > {log} 2>&1"


# ── Step 01: Compound collection ─────────────────────────────
rule step01_compounds:
    output: "data/processed/01_compounds_filtered.csv"
    log:    "logs/01_compounds.log"
    shell:
        "python scripts/python/01_compounds.py --config config.yaml > {log} 2>&1"


# ── Step 02: Target prediction ───────────────────────────────
rule step02_targets:
    input:  "data/processed/01_compounds_filtered.csv"
    output: "data/processed/02_targets_merged.csv"
    log:    "logs/02_targets.log"
    shell:
        "python scripts/python/02_targets.py --config config.yaml > {log} 2>&1"


# ── Step 03: Disease targets ──────────────────────────────────
rule step03_disease:
    output: "data/processed/03_disease_targets.csv"
    log:    "logs/03_disease.log"
    shell:
        "python scripts/python/03_disease.py --config config.yaml > {log} 2>&1"


# ── Step 04: Venn + PPI ───────────────────────────────────────
rule step04_ppi:
    input:
        drug    = "data/processed/02_targets_merged.csv",
        disease = "data/processed/03_disease_targets.csv"
    output:
        genes = "data/processed/04_intersection_genes.csv",
        edges = "data/processed/04_ppi_edges.csv",
        venn  = "results/figures/04_venn.pdf"
    log: "logs/04_ppi.log"
    shell:
        "python scripts/python/04_ppi.py --config config.yaml > {log} 2>&1"


# ── Step 05: Hub gene selection ───────────────────────────────
rule step05_hub_genes:
    input:
        "data/processed/04_intersection_genes.csv",
        "data/processed/04_ppi_edges.csv"
    output: "data/processed/05_hub_genes.csv"
    log:    "logs/05_hub_genes.log"
    shell:
        "python scripts/python/05_hub_genes.py --config config.yaml > {log} 2>&1"


# ── Step 06: Enrichment analysis ─────────────────────────────
rule step06_enrichment:
    input:  "data/processed/05_hub_genes.csv"
    output: "results/figures/06_enrichment_kegg.pdf"
    log:    "logs/06_enrichment.log"
    shell:
        "python scripts/python/06_enrichment.py --config config.yaml > {log} 2>&1"


# ── Step 07: Mendelian Randomization (R) ─────────────────────
rule step07_mr:
    input:  "data/processed/05_hub_genes.csv"
    output:
        results = "data/processed/07_mr_results.csv",
        forest  = touch("results/figures/07_mr_forest.pdf")
    log: "logs/07_mr.log"
    shell:
        "Rscript scripts/R/07_mr.R config.yaml > {log} 2>&1 || true"


# ── Step 08: Molecular docking ────────────────────────────────
rule step08_docking:
    input:
        compounds = "data/processed/01_compounds_filtered.csv",
        hub_genes = "data/processed/05_hub_genes.csv"
    output: "data/processed/08_docking_scores.csv"
    log:    "logs/08_docking.log"
    shell:
        "python scripts/python/08_docking.py --config config.yaml > {log} 2>&1"


# ── Step 09: scRNA-seq (optional) ────────────────────────────
rule step09_scrna:
    input:  "data/processed/05_hub_genes.csv"
    output: touch("data/processed/09_scrna_done.flag")
    log:    "logs/09_scrna.log"
    shell:
        "Rscript scripts/R/09_scrna.R config.yaml > {log} 2>&1"


rule step09b_immune:
    input:
        expression = "data/raw/cached/geo_expression_matrix.csv",
        hub_genes = "data/processed/05_hub_genes.csv"
    output:
        scores = "data/processed/09b_immune_scores.csv",
        heatmap = "results/figures/09b_immune_heatmap.pdf"
    log: "logs/09b_immune.log"
    shell:
        "python scripts/python/09b_immune.py --config config.yaml > {log} 2>&1"


# ── Step 10: Summary visualization ───────────────────────────
rule step10_visualization:
    input:
        enrichment = "results/figures/06_enrichment_kegg.pdf",
        mr         = "data/processed/07_mr_results.csv",
        docking    = "data/processed/08_docking_scores.csv"
    output: "results/figures/10_summary_figure.pdf"
    log:    "logs/10_visualization.log"
    shell:
        "python scripts/python/10_visualization.py --config config.yaml > {log} 2>&1"


# ── Utility rules ─────────────────────────────────────────────
rule clean:
    shell:
        "python -c \"from pathlib import Path; import shutil; "
        "targets=['data/processed','results/figures','results/tables','results/docking']; "
        "[shutil.rmtree(t, ignore_errors=True) for t in targets]; "
        "[p.unlink() for p in Path('logs').glob('*.log')]\""

rule clean_all:
    shell:
        "python -c \"from pathlib import Path; import shutil; "
        "targets=['data/processed','results','data/raw/cached']; "
        "[shutil.rmtree(t, ignore_errors=True) for t in ['data/processed','results']]; "
        "[p.unlink() for p in Path('data/raw/cached').glob('_*')] if Path('data/raw/cached').exists() else None; "
        "[p.unlink() for p in Path('logs').glob('*.log')]\""

rule dag:
    shell: "snakemake --dag | dot -Tpdf > docs/workflow_dag.pdf"

rule preflight:
    shell: "python scripts/python/00_preflight.py --config config.yaml"
