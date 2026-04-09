#!/usr/bin/env python3
"""
Step 05: Hub gene selection.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from utils_geo import (
    download_geo_matrix,
    extract_sample_labels,
    map_probes_to_genes,
    parse_geo_series_matrix,
)
from utils_validate import log_data_provenance, print_step_summary, validate_csv

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def load_config(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def compute_topology(edges_df: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """Compute degree and betweenness centrality from a PPI edge list."""
    import networkx as nx

    graph = nx.Graph()
    graph.add_nodes_from(genes)
    for _, row in edges_df.iterrows():
        gene1 = row["gene1"]
        gene2 = row["gene2"]
        if gene1 in genes and gene2 in genes:
            graph.add_edge(gene1, gene2, weight=row.get("combined_score", 1.0))

    degree = dict(graph.degree())
    betweenness = nx.betweenness_centrality(graph)
    topo = pd.DataFrame(
        {
            "gene_symbol": list(degree.keys()),
            "degree": list(degree.values()),
            "betweenness": [betweenness.get(gene, 0.0) for gene in degree],
        }
    )
    return topo.sort_values(["degree", "betweenness"], ascending=[False, False]).reset_index(drop=True)


def log_mcode_instructions(cfg: dict[str, Any]) -> None:
    mcode_cfg = cfg["hub_selection"]["mcode"]
    log.info("MCODE manual instructions:")
    log.info("  1. Import data/processed/04_ppi_network.sif into Cytoscape")
    log.info("  2. Run Apps -> MCODE -> Analyze current network")
    log.info("  3. Use degree_cutoff=%s", mcode_cfg["degree_cutoff"])
    log.info("  4. Use score_cutoff=%s", mcode_cfg["score_cutoff"])
    log.info("  5. Use k_core=%s", mcode_cfg["k_core"])
    log.info("  6. Save the node list to data/raw/cached/mcode_clusters.csv")


def load_mcode_genes(cache_dir: Path) -> set[str]:
    path = cache_dir / "mcode_clusters.csv"
    if not path.exists():
        return set()
    if not validate_csv(path, min_rows=1, step_name="Step 05 MCODE"):
        return set()

    df = pd.read_csv(path)
    gene_col = None
    for column in df.columns:
        lowered = column.lower()
        if any(token in lowered for token in ["gene", "symbol", "node", "name"]):
            gene_col = column
            break
    if not gene_col:
        log.warning("MCODE export did not contain a recognizable gene column.")
        return set()
    genes = set(df[gene_col].dropna().astype(str).str.upper().str.strip())
    log.info("Loaded %d MCODE genes.", len(genes))
    return genes


def _materialize_geo_expression(cache_dir: Path, cfg: dict[str, Any]) -> Path | None:
    expr_csv_path = cache_dir / "geo_expression_matrix.csv"
    if expr_csv_path.exists():
        return expr_csv_path

    accession = cfg["disease"].get("geo_bulk_accession")
    if not accession:
        log.info("No GEO bulk accession configured; skipping GEO expression bootstrap.")
        return None

    matrix_path = download_geo_matrix(accession, cache_dir)
    if matrix_path is None:
        return None

    try:
        expr_df, metadata_df = parse_geo_series_matrix(matrix_path)
        gene_expr = map_probes_to_genes(expr_df, metadata_df.attrs.get("platform_id", ""))
        if gene_expr.empty:
            log.warning("GEO expression mapping produced no gene-level matrix.")
            return None

        sample_expr = gene_expr.transpose()
        sample_expr.index.name = "sample_id"
        sample_expr.to_csv(expr_csv_path)
        log.info("Materialized GEO expression matrix -> %s", expr_csv_path)

        labels_df = extract_sample_labels(metadata_df)
        if labels_df is not None:
            labels_path = cache_dir / "geo_sample_labels.csv"
            labels_df.to_csv(labels_path, index=False)
            log.info("Materialized GEO sample labels -> %s", labels_path)
        else:
            log.warning("Could not infer case/control labels from GEO metadata automatically.")
        return expr_csv_path
    except Exception as exc:
        log.warning("Failed to parse GEO series matrix for %s: %s", accession, exc)
        return None


def load_geo_matrix(cache_dir: Path, genes: list[str], cfg: dict[str, Any]) -> pd.DataFrame | None:
    """
    Load or auto-materialize a GEO expression matrix.

    Expected final format:
      - rows = samples
      - cols = genes
    """
    expr_path = _materialize_geo_expression(cache_dir, cfg)
    if expr_path is None or not expr_path.exists():
        log.info("No GEO expression matrix available; using the network-only path.")
        return None

    df = pd.read_csv(expr_path, index_col=0)
    available = [gene for gene in genes if gene in df.columns]
    if len(available) < 5:
        log.warning("Only %d hub candidates overlap the GEO matrix; skipping ML layer.", len(available))
        return None
    log.info("Loaded GEO matrix: %d samples x %d genes", df.shape[0], len(available))
    return df[available]


def run_lasso(X: pd.DataFrame, y: np.ndarray, cfg: dict[str, Any]) -> pd.Series:
    from sklearn.linear_model import LassoCV
    from sklearn.preprocessing import StandardScaler

    lasso_cfg = cfg["hub_selection"]["lasso"]
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    alphas = np.logspace(
        np.log10(lasso_cfg["alpha_range"][0]),
        np.log10(lasso_cfg["alpha_range"][1]),
        50,
    )
    model = LassoCV(
        alphas=alphas,
        cv=lasso_cfg["cv_folds"],
        random_state=lasso_cfg["seed"],
        max_iter=10000,
    )
    model.fit(X_scaled, y)
    return pd.Series(model.coef_, index=X.columns, name="lasso_coef")


def run_random_forest(X: pd.DataFrame, y: np.ndarray, cfg: dict[str, Any]) -> pd.Series:
    from sklearn.ensemble import RandomForestRegressor

    rf_cfg = cfg["hub_selection"]["random_forest"]
    model = RandomForestRegressor(
        n_estimators=rf_cfg["n_estimators"],
        max_features=rf_cfg["max_features"],
        random_state=rf_cfg["seed"],
        n_jobs=-1,
    )
    model.fit(X, y)
    return pd.Series(model.feature_importances_, index=X.columns, name="rf_importance")


def _encode_labels(values: pd.Series | np.ndarray) -> np.ndarray:
    series = pd.Series(values)
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.notna().all():
        return numeric.to_numpy()
    categorical = pd.Categorical(series.astype(str))
    mapping = {category: code for code, category in enumerate(categorical.categories)}
    log.info("Encoded categorical phenotype labels: %s", mapping)
    return categorical.codes.astype(float)


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 05: Hub gene selection")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    hub_cfg = cfg["hub_selection"]

    processed_dir = Path("data/processed")
    cache_dir = Path("data/raw/cached")
    genes_path = processed_dir / "04_intersection_genes.csv"
    edges_path = processed_dir / "04_ppi_edges.csv"

    if not validate_csv(
        genes_path,
        required_columns=["gene_symbol"],
        min_rows=1,
        step_name="Step 05",
    ):
        sys.exit(1)
    if not validate_csv(
        edges_path,
        required_columns=["gene1", "gene2", "combined_score"],
        min_rows=1,
        step_name="Step 05",
    ):
        sys.exit(1)

    genes_df = pd.read_csv(genes_path)
    edges_df = pd.read_csv(edges_path)
    genes = genes_df["gene_symbol"].dropna().astype(str).str.upper().tolist()

    log.info("=== Step 05: Hub gene selection | %d intersection genes ===", len(genes))

    topo = compute_topology(edges_df, genes)
    top_degree = set(topo.nlargest(hub_cfg["top_n_degree"], "degree")["gene_symbol"])
    top_betweenness = set(topo.nlargest(hub_cfg["top_n_betweenness"], "betweenness")["gene_symbol"])
    topo_hub = top_degree | top_betweenness
    log.info("Topology layer retained %d candidate genes.", len(topo_hub))

    log_mcode_instructions(cfg)
    mcode_genes = load_mcode_genes(cache_dir)

    geo_matrix = load_geo_matrix(cache_dir, genes, cfg)
    ml_hub: set[str] = set()
    lasso_coef = pd.Series(dtype=float, name="lasso_coef")
    rf_importance = pd.Series(dtype=float, name="rf_importance")

    if geo_matrix is not None and geo_matrix.shape[0] >= 10 and geo_matrix.shape[1] >= 3:
        X = geo_matrix.copy()
        labels_path = cache_dir / "geo_sample_labels.csv"
        if labels_path.exists() and validate_csv(labels_path, min_rows=2, step_name="Step 05 GEO labels"):
            labels_df = pd.read_csv(labels_path)
            sample_col = labels_df.columns[0]
            label_col = labels_df.columns[1]
            labels_df = labels_df.rename(columns={sample_col: "sample_id", label_col: "label"}).set_index("sample_id")
            common_samples = X.index.intersection(labels_df.index)
            X = X.loc[common_samples]
            y = _encode_labels(labels_df.loc[common_samples, "label"])
            log.info("Using %d labeled GEO samples for the ML layer.", len(common_samples))
        else:
            from sklearn.decomposition import PCA

            y = PCA(n_components=1, random_state=42).fit_transform(X).ravel()
            log.info("No sample labels available; using expression PC1 as the ML target.")

        lasso_coef = run_lasso(X, y, cfg)
        rf_importance = run_random_forest(X, y, cfg)

        lasso_genes = set(lasso_coef[lasso_coef != 0].index)
        adaptive_threshold = max(
            float(hub_cfg["random_forest"]["importance_threshold"]),
            1.0 / max(len(rf_importance), 1),
        )
        rf_genes = set(rf_importance[rf_importance >= adaptive_threshold].index)
        ml_hub = lasso_genes & rf_genes
        log.info(
            "ML layer retained %d shared genes (LASSO=%d, RF=%d).",
            len(ml_hub),
            len(lasso_genes),
            len(rf_genes),
        )
    else:
        log.warning("ML layer skipped: no suitable GEO expression matrix was available.")

    vote: dict[str, int] = {}
    for gene in genes:
        vote[gene] = int(gene in topo_hub) + int(gene in mcode_genes) + int(gene in ml_hub)

    layers_active = 1 + int(bool(mcode_genes)) + int(bool(ml_hub))
    min_votes = min(2, layers_active)
    final_hub = {gene for gene, count in vote.items() if count >= min_votes}
    if len(final_hub) < 3:
        log.warning("Only %d genes passed the vote threshold; falling back to topology candidates.", len(final_hub))
        final_hub = topo_hub

    result = topo[topo["gene_symbol"].isin(final_hub)].copy()
    result = result.merge(
        lasso_coef.reset_index().rename(columns={"index": "gene_symbol"}),
        on="gene_symbol",
        how="left",
    )
    result = result.merge(
        rf_importance.reset_index().rename(columns={"index": "gene_symbol"}),
        on="gene_symbol",
        how="left",
    )
    result = result[["gene_symbol", "degree", "betweenness", "lasso_coef", "rf_importance"]]
    result = result.sort_values(["degree", "betweenness"], ascending=[False, False]).reset_index(drop=True)

    out_path = processed_dir / "05_hub_genes.csv"
    result.to_csv(out_path, index=False)
    log.info("Output -> %s (%d hub genes)", out_path, len(result))

    sources_used = ["topology"]
    if mcode_genes:
        sources_used.append("MCODE")
    if ml_hub:
        sources_used.append("GEO+ML")

    log_data_provenance(
        step="Step 05: Hub gene selection",
        output_path=out_path,
        sources=sources_used,
        cfg=cfg,
        extra={
            "n_input": int(len(genes)),
            "n_output": int(len(result)),
            "filters_applied": {
                "top_n_degree": hub_cfg["top_n_degree"],
                "top_n_betweenness": hub_cfg["top_n_betweenness"],
            },
        },
    )
    print_step_summary(
        "Step 05",
        out_path,
        cfg,
        next_step="python scripts/python/06_enrichment.py --config config.yaml",
    )


if __name__ == "__main__":
    main()
