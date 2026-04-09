#!/usr/bin/env python3
"""
Step 00: Preflight Check — Pipeline Readiness Validator
========================================================
Validates environment, dependencies, config, and data files
BEFORE running the pipeline. Reports clear pass/fail per item.

Usage:
  python scripts/python/00_preflight.py --config config.yaml
  python scripts/python/00_preflight.py --config config.yaml --step 02
"""

import argparse
import importlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

# ── ANSI colors ──────────────────────────────────────────────
GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
RESET = "\033[0m"
PASS = f"{GREEN}[PASS]{RESET}"
FAIL = f"{RED}[FAIL]{RESET}"
WARN = f"{YELLOW}[WARN]{RESET}"


def load_config(path: str) -> dict:
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


# ── Python dependency checks ────────────────────────────────
PYTHON_REQUIRED = [
    "pandas", "numpy", "requests", "yaml", "networkx",
    "sklearn", "matplotlib", "matplotlib_venn", "mygene",
]
PYTHON_OPTIONAL = ["rdkit", "openbabel"]


def check_python_deps():
    print("\n── Python Dependencies ──")
    all_ok = True
    for pkg in PYTHON_REQUIRED:
        try:
            importlib.import_module(pkg)
            print(f"  {PASS} {pkg}")
        except ImportError:
            print(f"  {FAIL} {pkg} — pip install {pkg}")
            all_ok = False
    for pkg in PYTHON_OPTIONAL:
        try:
            importlib.import_module(pkg)
            print(f"  {PASS} {pkg} (optional)")
        except ImportError:
            print(f"  {WARN} {pkg} (optional, needed for docking)")
    return all_ok


# ── R dependency checks ─────────────────────────────────────
R_REQUIRED = [
    "TwoSampleMR", "clusterProfiler", "org.Hs.eg.db", "enrichplot",
    "DOSE", "ggplot2", "dplyr", "readr", "yaml",
]
R_OPTIONAL = ["Seurat", "harmony", "GSVA", "Nebulosa", "survival", "survminer"]


def check_r_deps():
    print("\n── R Dependencies ──")
    rscript = shutil.which("Rscript")
    if not rscript:
        print(f"  {FAIL} Rscript not found in PATH")
        print(f"        Install R 4.3+ and add to PATH")
        return False

    print(f"  {PASS} Rscript found: {rscript}")
    all_ok = True

    for pkg in R_REQUIRED:
        cmd = f'cat(requireNamespace("{pkg}", quietly=TRUE))'
        try:
            r = subprocess.run(["Rscript", "-e", cmd], capture_output=True, text=True, timeout=15)
            if "TRUE" in r.stdout:
                print(f"  {PASS} {pkg}")
            else:
                print(f"  {FAIL} {pkg} — install in R: BiocManager::install('{pkg}')")
                all_ok = False
        except Exception:
            print(f"  {FAIL} {pkg} — check failed")
            all_ok = False

    for pkg in R_OPTIONAL:
        cmd = f'cat(requireNamespace("{pkg}", quietly=TRUE))'
        try:
            r = subprocess.run(["Rscript", "-e", cmd], capture_output=True, text=True, timeout=15)
            if "TRUE" in r.stdout:
                print(f"  {PASS} {pkg} (optional)")
            else:
                print(f"  {WARN} {pkg} (optional)")
        except Exception:
            print(f"  {WARN} {pkg} (optional, check failed)")

    return all_ok


# ── External tool checks ────────────────────────────────────
def check_external_tools():
    print("\n── External Tools ──")
    tools = {
        "vina": "AutoDock Vina (docking)",
        "obabel": "Open Babel (ligand preparation)",
        "gmx": "GROMACS (MD simulation)",
    }
    all_ok = True
    for cmd, desc in tools.items():
        if shutil.which(cmd):
            print(f"  {PASS} {desc}")
        else:
            print(f"  {WARN} {desc} — not in PATH (needed for Step 08)")
    return all_ok


# ── Config validation ────────────────────────────────────────
def check_config(cfg: dict):
    print("\n── Config Validation ──")
    ok = True
    markers = cfg.get("validation_rules", {}).get("placeholder_markers", [])

    # Check drug settings
    drug = cfg.get("drug", {})
    if not drug.get("name"):
        print(f"  {FAIL} drug.name is empty")
        ok = False
    else:
        print(f"  {PASS} drug.name = {drug['name']}")

    # Check disease settings
    disease = cfg.get("disease", {})
    if not disease.get("name"):
        print(f"  {FAIL} disease.name is empty")
        ok = False
    else:
        print(f"  {PASS} disease.name = {disease['name']}")

    if not disease.get("mesh_term"):
        print(f"  {WARN} disease.mesh_term not set (needed for DisGeNET)")
    else:
        print(f"  {PASS} disease.mesh_term = {disease['mesh_term']}")

    # Check GEO accessions
    geo_bulk = disease.get("geo_bulk_accession")
    geo_scrna = disease.get("geo_scrna_accession")
    if geo_bulk and not geo_bulk.startswith("GSE"):
        print(f"  {WARN} geo_bulk_accession '{geo_bulk}' doesn't look like a GEO ID")
    elif geo_bulk:
        print(f"  {PASS} geo_bulk_accession = {geo_bulk}")

    # Check query_date
    qd = cfg.get("target_prediction", {}).get("query_date", "")
    if qd in markers or "YYYY" in str(qd):
        print(f"  {WARN} target_prediction.query_date not filled (will use today's date)")
    else:
        print(f"  {PASS} query_date = {qd}")

    return ok


# ── Data file checks (per step) ─────────────────────────────
def check_data_files(cfg: dict, step: int | None = None):
    print("\n── Data Files ──")
    proc = Path("data/processed")
    cache = Path("data/raw/cached")
    markers = cfg.get("validation_rules", {}).get("placeholder_markers", [])
    rules = cfg.get("validation_rules", {})

    checks = [
        (1, proc / "01_compounds_filtered.csv", "gene_symbol|smiles|pubchem_cid|Molecule Name", rules.get("min_compounds", 5)),
        (2, proc / "02_targets_merged.csv", "gene_symbol", rules.get("min_drug_targets", 10)),
        (3, proc / "03_disease_targets.csv", "gene_symbol", rules.get("min_disease_targets", 50)),
        (4, proc / "04_intersection_genes.csv", "gene_symbol", rules.get("min_intersection", 5)),
        (4, proc / "04_ppi_edges.csv", "gene1", 1),
        (5, proc / "05_hub_genes.csv", "gene_symbol", rules.get("min_hub_genes", 3)),
        (7, proc / "07_mr_results.csv", "gene|method", 1),
        (8, proc / "08_docking_scores.csv", "compound|target", 1),
    ]

    for step_n, fpath, expected_cols, min_rows in checks:
        if step is not None and step_n != step:
            continue
        if not fpath.exists():
            print(f"  {WARN} Step {step_n:02d}: {fpath.name} — not yet generated")
            continue

        import pandas as pd
        try:
            df = pd.read_csv(fpath)
            n = len(df)

            # Check for placeholder content
            has_placeholder = False
            for col in df.columns:
                if df[col].dtype == object:
                    for marker in markers:
                        if df[col].str.contains(marker, case=False, na=False).any():
                            has_placeholder = True
                            break

            if has_placeholder:
                print(f"  {WARN} Step {step_n:02d}: {fpath.name} — {n} rows, CONTAINS PLACEHOLDERS")
            elif n < min_rows:
                print(f"  {FAIL} Step {step_n:02d}: {fpath.name} — {n} rows (min: {min_rows})")
            else:
                print(f"  {PASS} Step {step_n:02d}: {fpath.name} — {n} rows")
        except Exception as e:
            print(f"  {FAIL} Step {step_n:02d}: {fpath.name} — read error: {e}")

    # Check manual input files
    print("\n── Manual Input Files (data/raw/cached/) ──")
    manual_files = {
        "tcmsp_targets_manual.csv": "Step 02 — TCMSP herb-target pairs",
        "pharmmapper_manual.csv": "Step 02 — PharmMapper results",
        "ctd_manual.csv": "Step 02 — CTD interactions",
        "sea_manual.csv": "Step 02 — SEA results",
        "genecards_manual.csv": "Step 03 — GeneCards disease genes",
        "omim_manual.csv": "Step 03 — OMIM disease genes",
        "disgenet_manual.csv": "Step 03 — DisGeNET disease genes",
        "mcode_clusters.csv": "Step 05 — MCODE cluster genes",
        "geo_expression_matrix.csv": "Step 05 — GEO expression matrix (ML)",
        "geo_sample_labels.csv": "Step 05 — GEO sample labels (ML)",
    }
    for fname, desc in manual_files.items():
        fpath = cache / fname
        if fpath.exists():
            import pandas as pd
            try:
                df = pd.read_csv(fpath)
                print(f"  {PASS} {fname} — {len(df)} rows | {desc}")
            except Exception:
                print(f"  {WARN} {fname} — exists but unreadable | {desc}")
        else:
            print(f"  {WARN} {fname} — not found | {desc}")

    # Check eQTL and GWAS files for MR
    print("\n── MR Data Files ──")
    eqtl_files = list(cache.glob("eqtlgen_cis_*.tsv"))
    if eqtl_files:
        print(f"  {PASS} eQTLGen file: {eqtl_files[0].name}")
    else:
        print(f"  {WARN} eQTLGen cis-eQTL file not found")
        print(f"        Download: https://eqtlgen.org/cis-eqtls.html")
        print(f"        Save to: data/raw/cached/eqtlgen_cis_YYYYMMDD.tsv")

    gwas_files = list(cache.glob("finngen_r10_*.gz"))
    if gwas_files:
        print(f"  {PASS} FinnGen GWAS: {gwas_files[0].name}")
    else:
        print(f"  {WARN} FinnGen GWAS file not found")
        print(f"        Download: https://www.finngen.fi/en/access_results")
        print(f"        Phenotype: {cfg.get('mr', {}).get('gwas_phenotype', 'N/A')}")
        print(f"        Save to: data/raw/cached/finngen_r10_RA_YYYYMMDD.gz")


# ── Summary ──────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Preflight check for NP pipeline")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--step", type=int, default=None, help="Check specific step only")
    parser.add_argument("--skip-r", action="store_true", help="Skip R dependency checks")
    args = parser.parse_args()

    print("=" * 60)
    print("  Network Pharmacology Pipeline — Preflight Check")
    print("=" * 60)

    cfg = load_config(args.config)
    drug = cfg.get("drug", {}).get("name", "?")
    disease = cfg.get("disease", {}).get("name", "?")
    print(f"  Drug:    {drug}")
    print(f"  Disease: {disease}")

    py_ok = check_python_deps()
    r_ok = True
    if not args.skip_r:
        r_ok = check_r_deps()
    check_external_tools()
    cfg_ok = check_config(cfg)
    check_data_files(cfg, args.step)

    print("\n" + "=" * 60)
    if py_ok and r_ok and cfg_ok:
        print(f"  {PASS} Preflight checks passed. Pipeline ready.")
    else:
        print(f"  {WARN} Some checks failed. Review above and fix before running.")
    print("=" * 60)


if __name__ == "__main__":
    main()
