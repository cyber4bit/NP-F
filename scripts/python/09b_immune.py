#!/usr/bin/env python3
"""
Step 09b: ssGSEA immune infiltration analysis.
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


def write_immune_r(
    cfg: dict[str, Any],
    expr_path: Path,
    hub_path: Path,
    out_csv: Path,
    out_pdf: Path,
    script_path: Path,
) -> None:
    immune_cfg = cfg["immune"]
    script_path.write_text(
        f"""options(warn=1)
suppressPackageStartupMessages({{
  library(GSVA)
  library(GSEABase)
  library(readr)
}})

write_placeholder <- function(message) {{
  placeholder <- data.frame(
    hub_gene = "placeholder",
    immune_cell_type = "placeholder",
    spearman_r = NA_real_,
    p_value = NA_real_,
    immune_correlated = FALSE,
    status = "placeholder",
    stringsAsFactors = FALSE
  )
  write.csv(placeholder, "{out_csv}", row.names = FALSE)
  pdf("{out_pdf}", width = 7, height = 5)
  plot.new()
  text(0.5, 0.55, message, cex = 1.1)
  text(0.5, 0.40, "See logs for details.", cex = 0.9)
  dev.off()
}}

expr_file <- "{expr_path.as_posix()}"
hub_file <- "{hub_path.as_posix()}"
gmt_file <- "{immune_cfg['gene_sets_gmt']}"

if (!file.exists(expr_file)) {{
  write_placeholder("Immune infiltration not executed: missing expression matrix.")
  quit(save = "no", status = 0)
}}
if (!file.exists(hub_file)) {{
  write_placeholder("Immune infiltration not executed: missing hub gene file.")
  quit(save = "no", status = 0)
}}
if (!file.exists(gmt_file)) {{
  write_placeholder("Immune infiltration not executed: missing MSigDB C7 GMT file.")
  quit(save = "no", status = 0)
}}

tryCatch({{
  expr <- read.csv(expr_file, row.names = 1, check.names = FALSE)
  hub_genes <- read_csv(hub_file, show_col_types = FALSE)$gene_symbol
  hub_genes <- intersect(hub_genes, colnames(expr))
  if (length(hub_genes) < 2) {{
    write_placeholder("Immune infiltration not executed: too few hub genes in expression matrix.")
    quit(save = "no", status = 0)
  }}

  expr_mat <- t(as.matrix(expr))
  mode(expr_mat) <- "numeric"
  gene_sets <- getGmt(gmt_file)
  if (length(gene_sets) == 0) {{
    write_placeholder("Immune infiltration not executed: empty GMT file.")
    quit(save = "no", status = 0)
  }}

  ssgsea_scores <- gsva(expr_mat, gene_sets, method = "ssgsea", kcdf = "Gaussian")
  max_sets <- min(nrow(ssgsea_scores), {int(immune_cfg['n_immune_cell_types'])})
  if (max_sets < 1) {{
    write_placeholder("Immune infiltration not executed: no immune gene sets scored.")
    quit(save = "no", status = 0)
  }}
  ssgsea_scores <- ssgsea_scores[seq_len(max_sets), , drop = FALSE]

  rows <- list()
  idx <- 1
  for (gene in hub_genes) {{
    gene_vec <- as.numeric(expr[[gene]])
    for (cell_type in rownames(ssgsea_scores)) {{
      score_vec <- as.numeric(ssgsea_scores[cell_type, ])
      test <- suppressWarnings(cor.test(gene_vec, score_vec, method = "spearman"))
      rows[[idx]] <- data.frame(
        hub_gene = gene,
        immune_cell_type = cell_type,
        spearman_r = unname(test$estimate),
        p_value = test$p.value,
        immune_correlated = abs(unname(test$estimate)) > {float(immune_cfg['correlation_abs_r_min'])} &
          test$p.value < {float(immune_cfg['p_cutoff'])},
        status = "ok",
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }}
  }}

  result_df <- do.call(rbind, rows)
  write.csv(result_df, "{out_csv}", row.names = FALSE)

  cor_matrix <- xtabs(spearman_r ~ hub_gene + immune_cell_type, data = result_df)
  pdf(
    "{out_pdf}",
    width = max(8, ncol(cor_matrix) * 0.35 + 4),
    height = max(6, nrow(cor_matrix) * 0.4 + 3)
  )
  heatmap(
    cor_matrix,
    Rowv = NA,
    Colv = NA,
    scale = "none",
    col = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(101),
    margins = c(10, 10),
    main = "Hub gene vs immune infiltration (Spearman r)"
  )
  dev.off()
}}, error = function(e) {{
  write_placeholder(paste("Immune infiltration failed:", e$message))
  quit(save = "no", status = 0)
}})
""",
        encoding="utf-8",
    )


def run_r(script_path: Path) -> bool:
    try:
        result = subprocess.run(
            ["Rscript", str(script_path)],
            capture_output=True,
            text=True,
            timeout=1800,
        )
        if result.stdout.strip():
            for line in result.stdout.strip().splitlines():
                log.info("[R/09b] %s", line)
        if result.returncode != 0:
            stderr = result.stderr.strip() if result.stderr else ""
            log.error("R/09b failed: %s", stderr[:500])
            return False
        return True
    except FileNotFoundError:
        log.error("Rscript not found.")
        return False
    except subprocess.TimeoutExpired:
        log.error("R/09b immune analysis timed out.")
        return False


def main() -> None:
    parser = argparse.ArgumentParser(description="Step 09b: Immune infiltration analysis")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    cache_dir = Path("data/raw/cached")
    processed_dir = Path("data/processed")
    figure_dir = Path("results/figures")
    cache_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    expr_path = cache_dir / "geo_expression_matrix.csv"
    hub_path = processed_dir / "05_hub_genes.csv"
    if not validate_csv(expr_path, min_rows=1, step_name="Step 09b"):
        log.warning("Bulk expression matrix is missing or unreadable; placeholder outputs will be written.")
    if not validate_csv(hub_path, required_columns=["gene_symbol"], min_rows=1, step_name="Step 09b"):
        sys.exit(1)

    out_csv = processed_dir / "09b_immune_scores.csv"
    out_pdf = figure_dir / "09b_immune_heatmap.pdf"
    script_path = cache_dir / "_09b_immune.R"
    write_immune_r(cfg, expr_path, hub_path, out_csv, out_pdf, script_path)
    run_r(script_path)

    if out_csv.exists():
        df = pd.read_csv(out_csv)
        log_data_provenance(
            step="Step 09b: Immune infiltration analysis",
            output_path=out_csv,
            sources=["GSVA", "MSigDB C7"],
            cfg=cfg,
            extra={
                "n_input": int(pd.read_csv(hub_path).shape[0]),
                "n_output": int(len(df)),
                "filters_applied": cfg["immune"],
            },
        )
        print_step_summary(
            "Step 09b",
            out_csv,
            cfg,
            next_step="python scripts/python/10_visualization.py --config config.yaml",
        )
    else:
        log.error("Expected immune output CSV was not created.")
        sys.exit(1)


if __name__ == "__main__":
    main()
