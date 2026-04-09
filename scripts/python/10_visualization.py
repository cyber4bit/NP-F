#!/usr/bin/env python3
"""
Step 10: Summary Visualization
================================
Input:  data/processed/04_intersection_genes.csv
        data/processed/04_ppi_edges.csv
        data/processed/05_hub_genes.csv
        data/processed/07_mr_results.csv
        data/processed/08_docking_scores.csv
        config.yaml
Output: results/figures/10_summary_figure.pdf
        results/figures/10_ppi_network.pdf
        results/figures/10_hub_heatmap.pdf
        results/figures/10_mr_forest.pdf
        results/figures/10_docking_bar.pdf

Panels:
  1. Venn diagram (final version)
  2. PPI network (networkx + matplotlib)
  3. Hub gene heatmap
  4. MR forest plot (from 07_mr_results.csv)
  5. Docking binding energy bar chart
  6. CETSA / SPR curve template (CTS Fig 2G-J)
  7. Kaplan-Meier template (SPI Fig 2B, via R survminer)

Reference:
  CTS: Fig 2G-J (SPR sensorgram, CETSA melting curve)
  SPI: Fig 2B (Kaplan-Meier survival)
  BBR: Fig 5 (GROMACS RMSD/RMSF)

Query date: YYYY-MM-DD
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def load_config(path: str) -> dict:
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


# ── Panel 1: Venn diagram ──────────────────────────────────────────────────

def plot_venn(cfg: dict, out_dir: Path, fig_dir: Path) -> plt.Figure | None:
    """Final Venn diagram (drug ∩ disease targets)."""
    drug_path = out_dir / "02_targets_merged.csv"
    disease_path = out_dir / "03_disease_targets.csv"
    if not drug_path.exists() or not disease_path.exists():
        log.warning("Skipping Venn: missing Step 02/03 outputs.")
        return None
    try:
        from matplotlib_venn import venn2
    except ImportError:
        log.warning("matplotlib-venn not installed. pip install matplotlib-venn")
        return None

    drug_genes = set(pd.read_csv(drug_path)["gene_symbol"].dropna().str.upper())
    disease_genes = set(pd.read_csv(disease_path)["gene_symbol"].dropna().str.upper())

    fig, ax = plt.subplots(figsize=(6, 5))
    venn2([drug_genes, disease_genes],
          set_labels=(cfg["drug"]["name"], cfg["disease"]["name"]), ax=ax)
    ax.set_title("Drug–Disease Target Intersection")
    fig.savefig(fig_dir / "10_venn.pdf", dpi=300, bbox_inches="tight")
    log.info("Panel 1: Venn → 10_venn.pdf")
    return fig


# ── Panel 2: PPI network ──────────────────────────────────────────────────

def plot_ppi_network(out_dir: Path, fig_dir: Path) -> plt.Figure | None:
    """PPI network visualization with hub gene highlighting."""
    edges_path = out_dir / "04_ppi_edges.csv"
    hub_path = out_dir / "05_hub_genes.csv"
    if not edges_path.exists():
        log.warning("Skipping PPI network: missing 04_ppi_edges.csv")
        return None

    import networkx as nx

    edges = pd.read_csv(edges_path)
    G = nx.Graph()
    for _, row in edges.iterrows():
        G.add_edge(row["gene1"], row["gene2"],
                    weight=row.get("combined_score", 0.5))

    hub_genes = set()
    if hub_path.exists():
        hub_genes = set(pd.read_csv(hub_path)["gene_symbol"].dropna().str.upper())

    fig, ax = plt.subplots(figsize=(10, 10))
    pos = nx.spring_layout(G, seed=42, k=1.5/np.sqrt(max(G.number_of_nodes(), 1)))

    # Non-hub nodes
    non_hub = [n for n in G.nodes() if n not in hub_genes]
    nx.draw_networkx_nodes(G, pos, nodelist=non_hub, node_size=50,
                           node_color="#cccccc", alpha=0.6, ax=ax)
    # Hub nodes
    hub_in_g = [n for n in G.nodes() if n in hub_genes]
    nx.draw_networkx_nodes(G, pos, nodelist=hub_in_g, node_size=300,
                           node_color="#e74c3c", alpha=0.9, ax=ax)
    nx.draw_networkx_edges(G, pos, alpha=0.15, ax=ax)
    nx.draw_networkx_labels(G, pos, labels={n: n for n in hub_in_g},
                            font_size=7, font_weight="bold", ax=ax)

    ax.set_title(f"PPI Network ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)")
    ax.axis("off")
    fig.savefig(fig_dir / "10_ppi_network.pdf", dpi=300, bbox_inches="tight")
    log.info("Panel 2: PPI network → 10_ppi_network.pdf")
    return fig


# ── Panel 3: Hub gene heatmap ─────────────────────────────────────────────

def plot_hub_heatmap(out_dir: Path, fig_dir: Path) -> plt.Figure | None:
    """Heatmap of hub gene metrics (degree, betweenness, LASSO, RF)."""
    hub_path = out_dir / "05_hub_genes.csv"
    if not hub_path.exists():
        return None

    df = pd.read_csv(hub_path).set_index("gene_symbol")
    cols = [c for c in ["degree", "betweenness", "lasso_coef", "rf_importance"]
            if c in df.columns]
    if not cols:
        return None

    data = df[cols].fillna(0)
    # Normalize each column to 0-1 for visualization
    data_norm = (data - data.min()) / (data.max() - data.min() + 1e-9)

    fig, ax = plt.subplots(figsize=(6, max(3, len(data) * 0.35)))
    im = ax.imshow(data_norm.values, cmap="YlOrRd", aspect="auto")
    ax.set_yticks(range(len(data_norm)))
    ax.set_yticklabels(data_norm.index, fontsize=8)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=45, ha="right", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.6, label="Normalized score")
    ax.set_title("Hub Gene Metrics")
    fig.savefig(fig_dir / "10_hub_heatmap.pdf", dpi=300, bbox_inches="tight")
    log.info("Panel 3: Heatmap → 10_hub_heatmap.pdf")
    return fig


# ── Panel 4: MR forest plot ──────────────────────────────────────────────

def plot_mr_forest(out_dir: Path, fig_dir: Path) -> plt.Figure | None:
    """Forest plot from MR results (BBR template)."""
    mr_path = out_dir / "07_mr_results.csv"
    if not mr_path.exists():
        log.warning("Skipping MR forest: missing 07_mr_results.csv")
        return None

    df = pd.read_csv(mr_path)
    gene_col = "gene" if "gene" in df.columns else "exposure" if "exposure" in df.columns else None
    required = {"method", "b", "se"}
    missing = required - set(df.columns)
    if gene_col is None or missing:
        if gene_col is None:
            missing = missing | {"gene"}
        log.warning(f"MR results missing columns {missing}")
        return None

    # Use IVW results if available
    ivw = df[df["method"].str.contains("IVW", case=False, na=False)]
    plot_df = ivw if len(ivw) > 0 else df.head(20)
    plot_df = plot_df[~plot_df["method"].astype(str).str.contains("placeholder", case=False, na=False)]
    if plot_df.empty:
        log.warning("Skipping MR forest: only placeholder MR rows were found")
        return None

    fig, ax = plt.subplots(figsize=(8, max(3, len(plot_df) * 0.4)))
    y_pos = range(len(plot_df))
    b = plot_df["b"].values
    ci_lo = b - 1.96 * plot_df["se"].values
    ci_hi = b + 1.96 * plot_df["se"].values

    ax.errorbar(b, y_pos, xerr=1.96 * plot_df["se"].values,
                fmt="D", color="#2c3e50", markersize=5, capsize=3)
    ax.axvline(x=0, color="grey", linestyle="--", alpha=0.5)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(plot_df[gene_col].values, fontsize=8)
    ax.set_xlabel("Causal estimate (β)")
    ax.set_title("Mendelian Randomization — IVW")
    ax.invert_yaxis()
    fig.savefig(fig_dir / "10_mr_forest.pdf", dpi=300, bbox_inches="tight")
    log.info("Panel 4: MR forest → 10_mr_forest.pdf")
    return fig


# ── Panel 5: Docking bar chart ───────────────────────────────────────────

def plot_docking_bar(out_dir: Path, fig_dir: Path, cfg: dict) -> plt.Figure | None:
    """Binding energy bar chart."""
    dock_path = out_dir / "08_docking_scores.csv"
    if not dock_path.exists():
        return None

    df = pd.read_csv(dock_path).dropna(subset=["binding_energy_kcal_mol"])
    if df.empty:
        return None

    df = df.sort_values("binding_energy_kcal_mol")
    labels = df["compound"].astype(str) + " → " + df["target"].astype(str)
    cutoff = cfg["docking"]["binding_energy_cutoff"]

    fig, ax = plt.subplots(figsize=(8, max(3, len(df) * 0.35)))
    colors = ["#e74c3c" if e <= cutoff else "#95a5a6"
              for e in df["binding_energy_kcal_mol"]]
    ax.barh(range(len(df)), df["binding_energy_kcal_mol"], color=colors)
    ax.axvline(x=cutoff, color="black", linestyle="--", label=f"Cutoff ({cutoff})")
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Binding energy (kcal/mol)")
    ax.set_title("Molecular Docking Results")
    ax.legend()
    fig.savefig(fig_dir / "10_docking_bar.pdf", dpi=300, bbox_inches="tight")
    log.info("Panel 5: Docking → 10_docking_bar.pdf")
    return fig


# ── Panel 6: CETSA / SPR template ────────────────────────────────────────

def plot_cetsa_spr_template(fig_dir: Path, cfg: dict) -> plt.Figure:
    """Template figures for CETSA melting curve and SPR sensorgram (CTS Fig 2G-J)."""
    val = cfg.get("validation", {})
    spr_doses = val.get("SPR", {}).get("gradient_doses", [0.01, 0.1, 1, 10, 100])
    cetsa_temps = val.get("CETSA", {}).get("temperature_gradient", [40, 45, 50, 55, 60, 65])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # SPR sensorgram template
    ax = axes[0]
    t = np.linspace(0, 300, 500)
    for i, dose in enumerate(spr_doses):
        # Simulated association-dissociation curve
        ka = 0.01 * dose
        ru = dose * 10 * (1 - np.exp(-ka * t))
        ru[t > 150] = ru[t <= 150][-1] * np.exp(-0.005 * (t[t > 150] - 150))
        ax.plot(t, ru, label=f"{dose} μM")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Response (RU)")
    ax.set_title("SPR Sensorgram (template)")
    ax.legend(fontsize=7, title="Dose")

    # CETSA melting curve template
    ax = axes[1]
    temps = np.array(cetsa_temps)
    # Simulated melting: DMSO vs compound-treated
    dmso = 1.0 / (1 + np.exp(0.3 * (temps - 52)))
    treated = 1.0 / (1 + np.exp(0.3 * (temps - 57)))  # shifted Tm
    ax.plot(temps, dmso, "o-", label="DMSO", color="#3498db")
    ax.plot(temps, treated, "s-", label="Compound", color="#e74c3c")
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Relative band intensity")
    ax.set_title("CETSA Melting Curve (template)")
    ax.legend()

    fig.suptitle("Experimental Validation Templates (fill with real data)", fontsize=11)
    fig.tight_layout()
    fig.savefig(fig_dir / "10_cetsa_spr_template.pdf", dpi=300, bbox_inches="tight")
    log.info("Panel 6: CETSA/SPR template → 10_cetsa_spr_template.pdf")
    return fig


# ── Panel 7: Kaplan-Meier template (R survminer) ─────────────────────────

def generate_km_r_script(fig_dir: Path, cache_dir: Path):
    """Generate R script for KM survival plot (SPI Fig 2B template)."""
    script = cache_dir / "_10_km_template.R"
    script.write_text(f'''\
# Kaplan-Meier survival template (SPI paper Fig 2B)
# Requires: survival, survminer, and expression + clinical data
suppressPackageStartupMessages({{
  library(survival); library(survminer)
}})
# Replace with real clinical + expression data:
# clin <- read.csv("data/raw/cached/clinical_data.csv")
# clin$group <- ifelse(clin$hub_gene_expr > median(clin$hub_gene_expr), "High", "Low")
# fit <- survfit(Surv(time, status) ~ group, data=clin)
# pdf("{fig_dir}/10_km_survival.pdf", width=7, height=6)
# ggsurvplot(fit, data=clin, pval=TRUE, risk.table=TRUE,
#            palette=c("#e74c3c","#3498db"), title="Hub Gene Survival")
# dev.off()
cat("KM template generated. Provide clinical data to produce plot.\\n")
''')
    log.info(f"Panel 7: KM R template → {script}")
    try:
        subprocess.run(["Rscript", str(script)], capture_output=True, timeout=30)
    except Exception:
        pass


# ── Composite summary figure ─────────────────────────────────────────────

def compose_summary(figs: list[plt.Figure | None], fig_dir: Path):
    """Combine available panels into a multi-page PDF."""
    from matplotlib.backends.backend_pdf import PdfPages

    out_path = fig_dir / "10_summary_figure.pdf"
    with PdfPages(out_path) as pdf:
        for fig in figs:
            if fig is not None:
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
    log.info(f"Summary figure → {out_path} ({sum(f is not None for f in figs)} pages)")


# ── main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Step 10: Summary visualization")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = Path("data/processed")
    fig_dir = Path("results/figures")
    cache_dir = Path("data/raw/cached")
    fig_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    log.info("=== Step 10: Summary Visualization ===")

    figs = []

    # 1. Venn
    figs.append(plot_venn(cfg, out_dir, fig_dir))
    # 2. PPI network
    figs.append(plot_ppi_network(out_dir, fig_dir))
    # 3. Hub gene heatmap
    figs.append(plot_hub_heatmap(out_dir, fig_dir))
    # 4. MR forest plot
    figs.append(plot_mr_forest(out_dir, fig_dir))
    # 5. Docking bar chart
    figs.append(plot_docking_bar(out_dir, fig_dir, cfg))
    # 6. CETSA / SPR template
    figs.append(plot_cetsa_spr_template(fig_dir, cfg))
    # 7. KM template (R)
    generate_km_r_script(fig_dir, cache_dir)

    # Compose multi-page PDF
    compose_summary(figs, fig_dir)

    log.info("=" * 60)
    n_ok = sum(1 for f in figs if f is not None)
    log.info(f"Panels generated: {n_ok} / 7")
    log.info("Individual PDFs saved in results/figures/10_*.pdf")
    log.info("=" * 60)
    log.info("Pipeline complete! Review results in results/figures/")


if __name__ == "__main__":
    main()
