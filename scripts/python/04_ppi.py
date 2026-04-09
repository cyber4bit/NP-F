#!/usr/bin/env python3
"""
Step 04: Venn intersection and STRING PPI construction.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib_venn
import pandas as pd
import requests
import yaml

from utils_validate import log_data_provenance, print_step_summary, validate_csv

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def load_config(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def string_api_ppi(genes: list[str], species: int = 9606, confidence: float = 0.4) -> pd.DataFrame:
    """Fetch STRING network edges."""
    if not genes:
        return pd.DataFrame(columns=["gene1", "gene2", "combined_score"])

    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "%0d".join(genes),
        "species": species,
        "required_score": int(confidence * 1000),
        "caller_identity": "np_pipeline_opensource",
    }
    try:
        response = requests.post(url, data=params, timeout=30)
        response.raise_for_status()
        payload = response.json()
        if not payload:
            log.warning("STRING returned an empty network.")
            return pd.DataFrame(columns=["gene1", "gene2", "combined_score"])
        edges = pd.DataFrame(
            [
                {
                    "gene1": item["preferredName_A"],
                    "gene2": item["preferredName_B"],
                    "combined_score": item["score"],
                }
                for item in payload
            ]
        )
        return edges.drop_duplicates().reset_index(drop=True)
    except Exception as exc:
        log.error("STRING API error: %s", exc)
        return pd.DataFrame(columns=["gene1", "gene2", "combined_score"])


def plot_venn(
    drug_genes: set[str],
    disease_genes: set[str],
    out_path: Path,
    drug_name: str,
    disease_name: str,
) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    matplotlib_venn.venn2(
        [drug_genes, disease_genes],
        set_labels=(f"{drug_name}\nTargets", f"{disease_name}\nTargets"),
        ax=ax,
    )
    ax.set_title(f"Target Intersection: {drug_name} vs {disease_name}", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("Venn diagram -> %s", out_path)


def export_cytoscape_format(edges_df: pd.DataFrame, nodes: list[str], out_dir: Path) -> None:
    sif_path = out_dir / "04_ppi_network.sif"
    node_path = out_dir / "04_ppi_nodes.csv"

    with open(sif_path, "w", encoding="utf-8") as handle:
        for _, row in edges_df.iterrows():
            handle.write(f"{row['gene1']}\tpp\t{row['gene2']}\n")

    pd.DataFrame({"gene": nodes, "node_type": "intersection"}).to_csv(node_path, index=False)
    log.info("Cytoscape exports -> %s , %s", sif_path, node_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 04: Venn + PPI network")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    ppi_cfg = cfg["ppi"]
    drug_name = cfg["drug"]["name"]
    disease_name = cfg["disease"]["name"]

    processed_dir = Path("data/processed")
    figure_dir = Path("results/figures")
    figure_dir.mkdir(parents=True, exist_ok=True)

    drug_targets_path = processed_dir / "02_targets_merged.csv"
    disease_targets_path = processed_dir / "03_disease_targets.csv"
    if not validate_csv(
        drug_targets_path,
        required_columns=["gene_symbol", "source"],
        min_rows=1,
        step_name="Step 04",
    ):
        sys.exit(1)
    if not validate_csv(
        disease_targets_path,
        required_columns=["gene_symbol", "source"],
        min_rows=1,
        step_name="Step 04",
    ):
        sys.exit(1)

    drug_df = pd.read_csv(drug_targets_path)
    disease_df = pd.read_csv(disease_targets_path)
    drug_genes = set(drug_df["gene_symbol"].dropna().astype(str).str.upper())
    disease_genes = set(disease_df["gene_symbol"].dropna().astype(str).str.upper())
    intersection = sorted(drug_genes & disease_genes)

    log.info("Drug targets: %d", len(drug_genes))
    log.info("Disease targets: %d", len(disease_genes))
    log.info("Intersection genes: %d", len(intersection))

    venn_path = figure_dir / "04_venn.pdf"
    plot_venn(drug_genes, disease_genes, venn_path, drug_name, disease_name)

    intersection_df = pd.DataFrame({"gene_symbol": intersection})
    intersection_path = processed_dir / "04_intersection_genes.csv"
    intersection_df.to_csv(intersection_path, index=False)
    log.info("Intersection genes -> %s", intersection_path)

    edges_df = string_api_ppi(
        genes=intersection,
        species=int(ppi_cfg["species_id"]),
        confidence=float(ppi_cfg["confidence_score"]),
    )
    edges_path = processed_dir / "04_ppi_edges.csv"
    edges_df.to_csv(edges_path, index=False)
    log.info("PPI edges -> %s (%d row(s))", edges_path, len(edges_df))

    export_cytoscape_format(edges_df, intersection, processed_dir)

    input_sources = sorted(
        set(drug_df.get("source", pd.Series(dtype=str)).dropna().astype(str))
        | set(disease_df.get("source", pd.Series(dtype=str)).dropna().astype(str))
    )
    log_data_provenance(
        step="Step 04: Target intersection",
        output_path=intersection_path,
        sources=input_sources,
        cfg=cfg,
        extra={
            "n_input": int(len(drug_df) + len(disease_df)),
            "n_output": int(len(intersection_df)),
            "filters_applied": {},
        },
    )
    log_data_provenance(
        step="Step 04: STRING PPI network",
        output_path=edges_path,
        sources=["STRING"],
        cfg=cfg,
        extra={
            "n_input": int(len(intersection_df)),
            "n_output": int(len(edges_df)),
            "filters_applied": {"confidence_score": ppi_cfg["confidence_score"]},
        },
    )
    print_step_summary(
        "Step 04",
        intersection_path,
        cfg,
        next_step="python scripts/python/05_hub_genes.py --config config.yaml",
    )

    if not intersection:
        log.warning("No shared genes were found. Downstream hub-gene selection will fail loudly.")


if __name__ == "__main__":
    main()
