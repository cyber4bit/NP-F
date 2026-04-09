#!/usr/bin/env python3
"""
Step 06: Enrichment analysis wrapper.
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from utils_validate import log_data_provenance, print_step_summary, validate_csv

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def load_config(path: str) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def write_enrichment_r(
    genes: list[str],
    cfg: dict[str, Any],
    out_dir: Path,
    fig_dir: Path,
    script_path: Path,
) -> None:
    ec = cfg["enrichment"]
    gene_vec = 'c("' + '","'.join(genes) + '")'
    script_path.write_text(
        f"""options(warn=1)
suppressPackageStartupMessages({{
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
}})

genes <- {gene_vec}
tryCatch({{
  eg <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  cat("Mapped", nrow(eg), "of", length(genes), "genes\\n")

  go_all <- enrichGO(
    gene=eg$ENTREZID,
    OrgDb=org.Hs.eg.db,
    ont="ALL",
    pAdjustMethod="BH",
    pvalueCutoff={ec["p_cutoff"]},
    qvalueCutoff={ec["q_cutoff"]},
    minGSSize={ec["min_gene_set_size"]},
    maxGSSize={ec["max_gene_set_size"]},
    readable=TRUE
  )
  if (!is.null(go_all) && nrow(as.data.frame(go_all)) > 0) {{
    write.csv(as.data.frame(go_all), "{out_dir}/06_enrichment_go.csv", row.names=FALSE)
    pdf("{fig_dir}/06_enrichment_go.pdf", width=12, height=8)
    print(barplot(go_all, showCategory=10, split="ONTOLOGY") +
      facet_grid(ONTOLOGY ~ ., scale="free"))
    dev.off()
  }}

  kegg <- enrichKEGG(
    gene=eg$ENTREZID,
    organism="hsa",
    pAdjustMethod="BH",
    pvalueCutoff={ec["p_cutoff"]},
    qvalueCutoff={ec["q_cutoff"]},
    minGSSize={ec["min_gene_set_size"]},
    maxGSSize={ec["max_gene_set_size"]}
  )
  if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {{
    write.csv(as.data.frame(kegg), "{out_dir}/06_enrichment_kegg.csv", row.names=FALSE)
    pdf("{fig_dir}/06_enrichment_kegg.pdf", width=10, height=8)
    print(dotplot(kegg, showCategory=20, title="KEGG Pathway"))
    dev.off()
  }} else {{
    pdf("{fig_dir}/06_enrichment_kegg.pdf", width=5, height=5)
    plot.new()
    text(0.5, 0.5, "No significant KEGG pathways", cex=1.2)
    dev.off()
  }}
}}, error=function(e) {{
  cat("Enrichment failed:", e$message, "\\n")
  quit(save="no", status=1)
}})
""",
        encoding="utf-8",
    )


def write_gsva_r(genes: list[str], cfg: dict[str, Any], out_dir: Path, script_path: Path) -> None:
    gc = cfg["enrichment"]["gsva"]
    gene_vec = 'c("' + '","'.join(genes) + '")'
    script_path.write_text(
        f"""options(warn=1)
suppressPackageStartupMessages({{
  library(GSVA)
}})

expr_file <- "data/raw/cached/geo_expression_matrix.csv"
if (!file.exists(expr_file)) {{
  cat("GSVA skipped: place expression matrix at", expr_file, "\\n")
  quit(save="no", status=0)
}}

tryCatch({{
  expr <- as.matrix(read.csv(expr_file, row.names=1, check.names=FALSE))
  genes <- {gene_vec}
  genes <- intersect(genes, colnames(expr))
  if (length(genes) < 3) {{
    cat("Too few genes overlap the expression matrix\\n")
    quit(save="no", status=0)
  }}
  gsva_res <- gsva(t(expr), list(hub=genes), method="{gc["method"]}", kcdf="{gc["kcdf"]}")
  scores <- gsva_res["hub", ]
  hi <- quantile(scores, 1 - {gc["expression_top_pct"]}/100)
  lo <- quantile(scores, {gc["expression_top_pct"]}/100)
  grp <- ifelse(scores >= hi, "High", ifelse(scores <= lo, "Low", "Mid"))
  write.csv(
    data.frame(sample=names(scores), gsva=as.numeric(scores), group=grp),
    "{out_dir}/06_gsva_scores.csv",
    row.names=FALSE
  )
}}, error=function(e) {{
  cat("GSVA failed:", e$message, "\\n")
  quit(save="no", status=1)
}})
""",
        encoding="utf-8",
    )


def run_r(script_path: Path, label: str) -> bool:
    try:
        result = subprocess.run(
            ["Rscript", str(script_path)],
            capture_output=True,
            text=True,
            timeout=900,
        )
        if result.stdout.strip():
            for line in result.stdout.strip().splitlines():
                log.info("[R/%s] %s", label, line)
        if result.returncode != 0:
            stderr = result.stderr.strip() if result.stderr else ""
            log.error("R/%s failed: %s", label, stderr[:500])
            return False
        return True
    except FileNotFoundError:
        log.error("Rscript not found. Install R and add it to PATH.")
        return False
    except subprocess.TimeoutExpired:
        log.error("R/%s timed out.", label)
        return False


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 06: Enrichment analysis")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--skip_gsva", action="store_true")
    args = parser.parse_args()

    cfg = load_config(args.config)
    processed_dir = Path("data/processed")
    figure_dir = Path("results/figures")
    cache_dir = Path("data/raw/cached")
    figure_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    hub_path = processed_dir / "05_hub_genes.csv"
    inter_path = processed_dir / "04_intersection_genes.csv"
    if not validate_csv(hub_path, required_columns=["gene_symbol"], min_rows=1, step_name="Step 06"):
        sys.exit(1)
    if not validate_csv(inter_path, required_columns=["gene_symbol"], min_rows=1, step_name="Step 06"):
        sys.exit(1)

    hub_genes = pd.read_csv(hub_path)["gene_symbol"].dropna().astype(str).tolist()
    all_genes = pd.read_csv(inter_path)["gene_symbol"].dropna().astype(str).tolist()
    log.info("=== Step 06: Enrichment | %d genes ===", len(all_genes))

    enrich_script = cache_dir / "_06_enrichment.R"
    write_enrichment_r(all_genes, cfg, processed_dir, figure_dir, enrich_script)
    run_r(enrich_script, "enrichment")

    if not args.skip_gsva:
        gsva_script = cache_dir / "_06_gsva.R"
        write_gsva_r(hub_genes, cfg, processed_dir, gsva_script)
        run_r(gsva_script, "gsva")

    for csv_name, label in [
        ("06_enrichment_go.csv", "GO enrichment"),
        ("06_enrichment_kegg.csv", "KEGG enrichment"),
    ]:
        csv_path = processed_dir / csv_name
        if csv_path.exists():
            df = pd.read_csv(csv_path)
            log_data_provenance(
                step=f"Step 06: {label}",
                output_path=csv_path,
                sources=["clusterProfiler"],
                cfg=cfg,
                extra={"n_input": int(len(all_genes)), "n_output": int(len(df)), "filters_applied": cfg["enrichment"]},
            )

    kegg_path = processed_dir / "06_enrichment_kegg.csv"
    if kegg_path.exists():
        print_step_summary(
            "Step 06",
            kegg_path,
            cfg,
            next_step="Rscript scripts/R/07_mr.R config.yaml",
        )
    else:
        log.warning("KEGG result CSV was not produced.")


if __name__ == "__main__":
    main()
