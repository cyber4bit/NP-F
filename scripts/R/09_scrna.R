#!/usr/bin/env Rscript
# ============================================================
# Step 09: scRNA-seq + Spatial Transcriptomics Analysis
# ============================================================
# Input:  data/raw/geo_datasets/<GSE_accession>/
#         data/processed/05_hub_genes.csv
#         config.yaml → scrna settings
# Output: results/figures/09_umap_clusters.pdf
#         results/figures/09_hub_target_dotplot.pdf
#         results/figures/09_spatial_featureplot.pdf  (if spatial)
#         data/processed/09_cell_type_expression.csv
#
# Method: Seurat 4.3 + Harmony batch correction (CTS template)
#         Nebulosa density estimation for sparse signal (SPI template)
#         ssGSEA immune infiltration (CTS: GSE141549, 28 cell types)
#         Spatial: Visium SpatialFeaturePlot + Spearman correlation
#
# Reference:
#   CTS paper: GSE213216 (scRNA), GSE141549 (bulk + ssGSEA)
#              Seurat 4.3, Harmony, resolution 0.1 (global), 0.2 (T/NK)
#              QC: min 500, max 4500 genes, < 5% mito
#   SPI paper: GSE141445 + GSE176031, Nebulosa KDE,
#              Visium (Sparkle DB), linkET Spearman, TIP platform
#
# Query date: YYYY-MM-DD
# GEO accession: see config.yaml
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(yaml)
  library(readr)
  library(patchwork)
})

# Optional: Nebulosa for sparse signal visualization (SPI template)
if (requireNamespace("Nebulosa", quietly = TRUE)) {
  library(Nebulosa)
  HAS_NEBULOSA <- TRUE
} else {
  message("Nebulosa not installed. Skipping KDE visualization.")
  HAS_NEBULOSA <- FALSE
}

# ── Config ───────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else "config.yaml"
cfg <- yaml::read_yaml(config_path)
sc_cfg <- cfg$scrna

if (!sc_cfg$enabled) {
  cat("scRNA-seq module disabled (scrna.enabled: false in config.yaml).\n")
  cat("Set scrna.enabled: true and provide GEO accession to run.\n")
  quit(status = 0)
}

set.seed(cfg$output$seed_global)
cat("=== Step 09: scRNA-seq / Spatial Analysis ===\n")
cat(sprintf("GEO scRNA accession:   %s\n", cfg$disease$geo_scrna_accession))
cat(sprintf("GEO bulk accession:    %s\n", cfg$disease$geo_bulk_accession))

# ── Load hub targets ──────────────────────────────────────────
hub_path <- "data/processed/05_hub_genes.csv"
if (!file.exists(hub_path)) stop("Missing: ", hub_path, ". Run Step 05 first.")
hub_genes <- read_csv(hub_path, show_col_types = FALSE)$gene_symbol
cat(sprintf("Hub genes to map:  %d\n", length(hub_genes)))

# ── Helper: load 10x Genomics matrix ─────────────────────────
load_10x_data <- function(data_dir, sample_name) {
  cat(sprintf("  Loading sample: %s\n", sample_name))
  counts <- Read10X(data.dir = data_dir)
  CreateSeuratObject(
    counts   = counts,
    project  = sample_name,
    min.cells    = 3,
    min.features = 200
  )
}

# ── Quality control ───────────────────────────────────────────
run_qc <- function(seurat_obj, cfg) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  # Filters (CTS template)
  seurat_obj <- subset(seurat_obj,
    subset = nFeature_RNA >= cfg$scrna$min_genes_per_cell &
             nFeature_RNA <= cfg$scrna$max_genes_per_cell &
             percent.mt   <  cfg$scrna$max_mito_pct
  )
  cat(sprintf("  After QC: %d cells\n", ncol(seurat_obj)))
  seurat_obj
}

# ── Standard Seurat pipeline ──────────────────────────────────
run_seurat_pipeline <- function(seurat_obj, dims = 20, resolution = 0.1) {
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:dims, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, verbose = FALSE)
  seurat_obj
}

# ── Harmony batch correction (CTS: multiple samples) ─────────
run_harmony <- function(seurat_list, dims = 20, resolution = 0.1) {
  cat("  Merging samples for Harmony integration...\n")
  merged <- merge(seurat_list[[1]], y = seurat_list[-1],
                  add.cell.ids = names(seurat_list))
  merged <- NormalizeData(merged, verbose = FALSE)
  merged <- FindVariableFeatures(merged, verbose = FALSE)
  merged <- ScaleData(merged, verbose = FALSE)
  merged <- RunPCA(merged, verbose = FALSE)

  # Harmony integration
  merged <- RunHarmony(merged, group.by.vars = "orig.ident",
                       max.iter.harmony = 20, verbose = FALSE)
  merged <- RunUMAP(merged, reduction = "harmony", dims = 1:dims, verbose = FALSE)
  merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:dims, verbose = FALSE)
  merged <- FindClusters(merged, resolution = resolution, verbose = FALSE)
  merged
}

# ── ssGSEA immune infiltration (CTS: bulk GSE141549) ─────────
run_ssgsea <- function(expr_matrix, gene_sets_path = NULL) {
  # Requires GSVA package
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    message("GSVA not installed. Skipping ssGSEA.")
    return(NULL)
  }
  library(GSVA)
  library(GSEABase)

  if (is.null(gene_sets_path)) {
    # Default: download immune cell gene sets from msigDB
    # Using GSVA built-in gene sets as placeholder
    message("Provide immune cell gene sets path for ssGSEA (MSigDB C7 or custom).")
    return(NULL)
  }

  gene_sets <- getGmt(gene_sets_path)
  gsva(expr_matrix, gene_sets, method = "ssgsea", verbose = FALSE)
}

# ── Main scRNA workflow ───────────────────────────────────────
geo_accession <- cfg$disease$geo_scrna_accession
geo_dir       <- file.path("data/raw/geo_datasets", geo_accession)

if (!dir.exists(geo_dir)) {
  cat(sprintf("\n⚠ GEO dataset not found: %s\n", geo_dir))
  cat("  To download:\n")
  cat(sprintf("    mkdir -p %s\n", geo_dir))
  cat(sprintf("    # Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s\n",
              geo_accession))
  cat("  Place 10x matrix files (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)\n")
  cat(sprintf("  or pre-processed RDS at: %s/seurat_object.rds\n", geo_dir))
  quit(status = 0)
}

# Check for pre-processed RDS (saves time on re-runs)
rds_path <- file.path(geo_dir, "seurat_object.rds")
if (file.exists(rds_path)) {
  cat(sprintf("Loading pre-processed Seurat object: %s\n", rds_path))
  seurat_merged <- readRDS(rds_path)
} else {
  # Load raw 10x data
  sample_dirs <- list.dirs(geo_dir, recursive = FALSE)
  if (length(sample_dirs) == 0) {
    # Single sample
    sample_dirs <- geo_dir
  }

  seurat_list <- lapply(seq_along(sample_dirs), function(i) {
    obj <- load_10x_data(sample_dirs[i], basename(sample_dirs[i]))
    run_qc(obj, cfg)
  })
  names(seurat_list) <- basename(sample_dirs)

  # Integration
  if (length(seurat_list) > 1) {
    seurat_merged <- run_harmony(seurat_list, dims = sc_cfg$harmony_dims,
                                 resolution = sc_cfg$umap_resolution_global)
  } else {
    seurat_merged <- run_seurat_pipeline(seurat_list[[1]],
                                         dims = sc_cfg$harmony_dims,
                                         resolution = sc_cfg$umap_resolution_global)
  }

  # Cache processed object
  saveRDS(seurat_merged, rds_path)
  cat(sprintf("Seurat object cached → %s\n", rds_path))
}

# ── UMAP + cell type annotation ───────────────────────────────
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

p_umap <- DimPlot(seurat_merged, reduction = "umap", label = TRUE, pt.size = 0.3) +
  ggtitle(sprintf("UMAP: %s | %s", cfg$disease$name, geo_accession)) +
  theme_bw()
ggsave("results/figures/09_umap_clusters.pdf", p_umap, width = 8, height = 7)
cat("UMAP → results/figures/09_umap_clusters.pdf\n")

# ── Hub target expression across clusters ─────────────────────
# Filter to genes present in the dataset
hub_present <- hub_genes[hub_genes %in% rownames(seurat_merged)]
cat(sprintf("Hub genes in dataset: %d / %d\n", length(hub_present), length(hub_genes)))

if (length(hub_present) > 0) {
  # Dot plot (CTS: Bubble plot of markers, Fig 1C/E)
  p_dot <- DotPlot(seurat_merged, features = hub_present) +
    RotatedAxis() +
    ggtitle("Hub Target Expression Across Cell Clusters") +
    theme_bw()
  ggsave("results/figures/09_hub_target_dotplot.pdf", p_dot,
         width = max(8, length(hub_present) * 0.6 + 2), height = 7)
  cat("Dot plot → results/figures/09_hub_target_dotplot.pdf\n")

  # Nebulosa KDE for sparse signal (SPI template)
  if (HAS_NEBULOSA && length(hub_present) >= 1) {
    top_target <- hub_present[1]  # visualize top hub gene
    tryCatch({
      p_kde <- plot_density(seurat_merged, features = top_target) +
        ggtitle(sprintf("Nebulosa KDE: %s", top_target))
      ggsave(sprintf("results/figures/09_nebulosa_%s.pdf", top_target), p_kde,
             width = 6, height = 5)
      cat(sprintf("Nebulosa KDE → results/figures/09_nebulosa_%s.pdf\n", top_target))
    }, error = function(e) message("Nebulosa failed: ", e$message))
  }

  # Cell type-resolved expression table
  avg_expr <- AverageExpression(seurat_merged, features = hub_present,
                                 group.by = "seurat_clusters", slot = "data")$RNA
  expr_df <- as.data.frame(t(avg_expr))
  expr_df$cluster <- rownames(expr_df)
  write_csv(expr_df, "data/processed/09_cell_type_expression.csv")
  cat("Cell type expression → data/processed/09_cell_type_expression.csv\n")
}

# ── Spatial transcriptomics (if enabled) ─────────────────────
if (isTRUE(sc_cfg$spatial$enabled)) {
  spatial_accession <- cfg$disease$geo_spatial_accession
  if (is.null(spatial_accession)) {
    message("No spatial GEO accession set (geo_spatial_accession: null). Skipping.")
  } else {
    spatial_dir <- file.path("data/raw/geo_datasets", spatial_accession)
    if (dir.exists(spatial_dir)) {
      cat(sprintf("Loading Visium spatial data: %s\n", spatial_accession))
      # Standard Seurat Visium loading
      visium <- Load10X_Spatial(spatial_dir)
      visium <- NormalizeData(visium)

      if (length(hub_present) > 0) {
        # SpatialFeaturePlot (SPI template: Fig 5B,F)
        p_spatial <- SpatialFeaturePlot(visium, features = hub_present[1]) +
          ggtitle(sprintf("Spatial: %s | %s", hub_present[1], spatial_accession))
        ggsave("results/figures/09_spatial_featureplot.pdf", p_spatial,
               width = 6, height = 5)
        cat("Spatial feature plot → results/figures/09_spatial_featureplot.pdf\n")
      }
    } else {
      message("Spatial data directory not found: ", spatial_dir)
    }
  }
}

cat("\n=== Step 09 complete ===\n")
cat("NEXT STEP: python scripts/python/10_visualization.py --config config.yaml\n")
