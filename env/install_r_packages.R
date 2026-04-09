#!/usr/bin/env Rscript
# ============================================================
# Install R packages for NP pipeline
# Usage: Rscript env/install_r_packages.R
# ============================================================

cat("Installing R packages for Network Pharmacology Pipeline...\n")

# ── CRAN packages ──
cran_pkgs <- c("yaml", "dplyr", "readr", "ggplot2", "patchwork",
               "survival", "survminer", "devtools")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s from CRAN...\n", pkg))
    install.packages(pkg, repos = "https://cran.r-project.org")
  } else {
    cat(sprintf("  OK: %s\n", pkg))
  }
}

# ── Bioconductor packages ──
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.r-project.org")
}

bioc_pkgs <- c(
  # Step 06: Enrichment
  "clusterProfiler", "org.Hs.eg.db", "DOSE", "enrichplot",
  # Step 06: GSVA
  "GSVA", "GSEABase",
  # Step 07: MR
  "TwoSampleMR",
  # Step 09: scRNA
  "Seurat", "harmony", "Nebulosa"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s from Bioconductor...\n", pkg))
    tryCatch(
      BiocManager::install(pkg, ask = FALSE, update = FALSE),
      error = function(e) cat(sprintf("  FAILED: %s — %s\n", pkg, e$message))
    )
  } else {
    cat(sprintf("  OK: %s\n", pkg))
  }
}

# ── TwoSampleMR (GitHub) ──
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  cat("Installing TwoSampleMR from GitHub...\n")
  tryCatch(
    devtools::install_github("MRCIEU/TwoSampleMR", upgrade = "never"),
    error = function(e) cat(sprintf("  FAILED: TwoSampleMR — %s\n", e$message))
  )
}

cat("\n=== R package installation complete ===\n")
cat("Run: Rscript -e 'installed.packages()[,c(1,3)]' to verify\n")
