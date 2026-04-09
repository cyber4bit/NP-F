#!/usr/bin/env Rscript
options(warn = 1)

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(yaml)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else "config.yaml"
cfg <- yaml::read_yaml(config_path)
mr_cfg <- cfg$mr

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

create_message_pdf <- function(path, message) {
  pdf(path, width = 7, height = 5)
  plot.new()
  text(0.5, 0.58, message, cex = 1.0)
  text(0.5, 0.40, "See logs/07_mr.log for details.", cex = 0.9)
  dev.off()
}

write_provenance <- function(n_input, n_output) {
  prov <- list(
    step = "Step 07: Mendelian Randomization",
    output_file = "07_mr_results.csv",
    generated_date = as.character(Sys.Date()),
    drug = cfg$drug$name,
    disease = cfg$disease$name,
    sources = c(cfg$mr$eqtl_source, cfg$mr$gwas_source),
    config_snapshot = list(seed = cfg$output$seed_global),
    n_input = n_input,
    n_output = n_output,
    filters_applied = list(
      snp_p_threshold = cfg$mr$snp_p_threshold,
      r2_clumping = cfg$mr$r2_clumping,
      f_stat_min = cfg$mr$f_stat_min
    )
  )
  writeLines(
    yaml::as.yaml(prov),
    "data/processed/07_mr_results.provenance.yaml"
  )
}

write_placeholder_outputs <- function(genes, message) {
  placeholder <- data.frame(
    gene = genes,
    method = "placeholder",
    b = NA_real_,
    se = NA_real_,
    pval = NA_real_,
    or = NA_real_,
    or_lci95 = NA_real_,
    or_uci95 = NA_real_,
    stringsAsFactors = FALSE
  )
  write_csv(placeholder, "data/processed/07_mr_results.csv")
  write_provenance(length(genes), nrow(placeholder))
  create_message_pdf("results/figures/07_mr_forest.pdf", message)
  cat("Placeholder -> data/processed/07_mr_results.csv\n")
  cat("Placeholder PDF -> results/figures/07_mr_forest.pdf\n")
}

safe_read_tsv <- function(path) {
  tryCatch(
    read_tsv(path, show_col_types = FALSE),
    error = function(e) {
      message("Failed to read ", path, ": ", e$message)
      NULL
    }
  )
}

load_eqtl_instruments <- function(genes, cfg) {
  eqtl_files <- Sys.glob("data/raw/cached/eqtlgen_cis_*.tsv")
  if (length(eqtl_files) == 0) {
    message("eQTLGen cis-eQTL file not found.")
    message("  Download from https://eqtlgen.org/cis-eqtls.html")
    return(NULL)
  }

  eqtl <- safe_read_tsv(eqtl_files[1])
  if (is.null(eqtl)) {
    return(NULL)
  }

  required_cols <- c("SNP", "beta", "se", "pval", "effect_allele",
                     "other_allele", "eaf", "gene_symbol", "N")
  missing_cols <- setdiff(required_cols, names(eqtl))
  if (length(missing_cols) > 0) {
    message("eQTL file is missing columns: ", paste(missing_cols, collapse = ", "))
    return(NULL)
  }

  eqtl %>%
    filter(gene_symbol %in% genes, pval < as.numeric(cfg$mr$snp_p_threshold)) %>%
    mutate(F_stat = (beta / se)^2) %>%
    filter(F_stat > cfg$mr$f_stat_min)
}

load_disease_gwas <- function() {
  gwas_files <- Sys.glob("data/raw/cached/finngen_r10_*.gz")
  if (length(gwas_files) == 0) {
    message("FinnGen GWAS file not found.")
    message("  Download from https://www.finngen.fi/en/access_results")
    return(NULL)
  }

  gwas <- safe_read_tsv(gwas_files[1])
  if (is.null(gwas)) {
    return(NULL)
  }

  required_cols <- c("rsid", "beta", "sebeta", "pval", "alt", "ref", "af_alt")
  missing_cols <- setdiff(required_cols, names(gwas))
  if (length(missing_cols) > 0) {
    message("GWAS file is missing columns: ", paste(missing_cols, collapse = ", "))
    return(NULL)
  }
  gwas
}

run_mr_analysis <- function(exposure_dat, outcome_dat, gene) {
  tryCatch({
    harmonized <- harmonise_data(exposure_dat, outcome_dat)
    if (nrow(harmonized) == 0) {
      message("No harmonized SNPs for ", gene)
      return(NULL)
    }

    methods <- if (nrow(harmonized) == 1) {
      c("mr_wald_ratio")
    } else {
      c(
        "mr_ivw",
        "mr_egger_regression",
        "mr_weighted_median",
        "mr_simple_mode",
        "mr_weighted_mode"
      )
    }

    res <- mr(harmonized, method_list = methods)
    res$gene <- gene
    res
  }, error = function(e) {
    message("MR failed for ", gene, ": ", e$message)
    NULL
  })
}

hub_path <- "data/processed/05_hub_genes.csv"
if (!file.exists(hub_path)) {
  stop("Missing data/processed/05_hub_genes.csv. Run Step 05 first.")
}

hub_genes <- tryCatch(
  read_csv(hub_path, show_col_types = FALSE),
  error = function(e) stop("Failed to read hub genes: ", e$message)
)
gene_list <- unique(na.omit(hub_genes$gene_symbol))

cat("=== Step 07: Mendelian Randomization ===\n")
cat(sprintf("Drug: %s\n", cfg$drug$name))
cat(sprintf("Disease: %s\n", cfg$disease$name))
cat(sprintf("Hub genes: %d\n", length(gene_list)))

iv_data <- load_eqtl_instruments(gene_list, cfg)
gwas_data <- load_disease_gwas()

if (is.null(iv_data) || is.null(gwas_data) || nrow(iv_data) == 0) {
  write_placeholder_outputs(
    gene_list,
    "MR not executed - missing input data. See logs/07_mr.log."
  )
  quit(save = "no", status = 0)
}

all_results <- list()
for (gene in gene_list) {
  gene_iv <- iv_data %>% filter(gene_symbol == gene)
  if (nrow(gene_iv) == 0) {
    message("No instruments found for ", gene)
    next
  }

  exposure_dat <- tryCatch(
    format_data(
      gene_iv,
      type = "exposure",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      eaf_col = "eaf",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      pval_col = "pval",
      samplesize_col = "N"
    ),
    error = function(e) {
      message("Failed to format exposure data for ", gene, ": ", e$message)
      NULL
    }
  )
  if (is.null(exposure_dat) || nrow(exposure_dat) == 0) {
    next
  }
  exposure_dat$exposure <- gene

  exposure_dat <- tryCatch(
    clump_data(exposure_dat, clump_r2 = mr_cfg$r2_clumping),
    error = function(e) {
      message("Clumping skipped for ", gene, ": ", e$message)
      exposure_dat
    }
  )

  outcome_dat <- tryCatch(
    format_data(
      gwas_data,
      type = "outcome",
      snp_col = "rsid",
      beta_col = "beta",
      se_col = "sebeta",
      pval_col = "pval",
      effect_allele_col = "alt",
      other_allele_col = "ref",
      eaf_col = "af_alt"
    ),
    error = function(e) {
      message("Failed to format outcome data: ", e$message)
      NULL
    }
  )
  if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
    next
  }
  outcome_dat$outcome <- cfg$disease$name

  res <- run_mr_analysis(exposure_dat, outcome_dat, gene)
  if (!is.null(res)) {
    all_results[[gene]] <- res
  }
}

if (length(all_results) == 0) {
  write_placeholder_outputs(
    gene_list,
    "MR not executed - no valid instruments remained after filtering."
  )
  quit(save = "no", status = 0)
}

combined <- bind_rows(all_results) %>%
  mutate(
    gene = ifelse(is.na(gene), exposure, gene),
    or = exp(b),
    or_lci95 = exp(b - 1.96 * se),
    or_uci95 = exp(b + 1.96 * se)
  ) %>%
  select(gene, method, b, se, pval, or, or_lci95, or_uci95)

write_csv(combined, "data/processed/07_mr_results.csv")
write_provenance(length(gene_list), nrow(combined))
cat(sprintf("MR results -> data/processed/07_mr_results.csv (%d rows)\n", nrow(combined)))

ivw_res <- combined %>%
  filter(grepl("inverse variance weighted|ivw", method, ignore.case = TRUE)) %>%
  arrange(pval)

if (nrow(ivw_res) > 0) {
  p_forest <- ggplot(
    ivw_res,
    aes(x = b, y = reorder(gene, -pval), xmin = b - 1.96 * se, xmax = b + 1.96 * se)
  ) +
    geom_point(size = 3, color = "#2c3e50") +
    geom_errorbarh(height = 0.3, color = "#2c3e50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = paste("MR Forest Plot:", cfg$drug$name, "to", cfg$disease$name),
      subtitle = "Method: Inverse variance weighted",
      x = "Beta (causal estimate)",
      y = "Gene"
    ) +
    theme_bw(base_size = 12)

  ggsave(
    "results/figures/07_mr_forest.pdf",
    p_forest,
    width = 8,
    height = max(4, nrow(ivw_res) * 0.4 + 2)
  )
  cat("Forest plot -> results/figures/07_mr_forest.pdf\n")
} else {
  create_message_pdf(
    "results/figures/07_mr_forest.pdf",
    "MR executed, but no IVW result was available for plotting."
  )
}

cat("=== Step 07 complete ===\n")
