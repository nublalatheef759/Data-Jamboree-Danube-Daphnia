## =============================================================================
## MOFA2 Multi-Omics Integration — Data Jamboree
## Role: Multi-omics integration (transcriptomics + metabolomics pos/neg)
## =============================================================================

# ── PART 0: Installation (run once, then comment out) ────────────────────────

# Install MOFA2 from Bioconductor

# Install Python backend (MOFA2 needs this)
# If reticulate asks you to install Miniconda, say YES
reticulate::install_miniconda()  # only if you don't have conda/miniconda
reticulate::py_install("mofapy2", pip = TRUE)

# Other packages

# ── PART 1: Load libraries ──────────────────────────────────────────────────

library(MOFA2)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# ── PART 2: Set your working directory ──────────────────────────────────────

# CHANGE THIS to where your files are

# ── PART 3: Load the three input matrices ───────────────────────────────────

# Transcriptomics: top 5000 most variable VST genes
rna <- read.csv("vst_top5000_for_iDEP.csv", row.names = 1, check.names = FALSE)
cat("RNA dimensions:", dim(rna), "\n")

# Metabolomics negative mode: PQN + imputed + glog
metabo_neg <- read.csv("polar_neg_pqn_imputed_glog.csv", row.names = 1, check.names = FALSE)
cat("Metabo Neg dimensions:", dim(metabo_neg), "\n")

# Metabolomics positive mode: PQN + imputed + glog  *** THE CORRECT FILE ***
# Haziq used polar_pos_pqn.csv (no glog) — that was the error
metabo_pos <- read.csv("polar_pos_pqn_imputed_glog.csv", row.names = 1, check.names = FALSE)
cat("Metabo Pos dimensions:", dim(metabo_pos), "\n")

# ── PART 4: Check orientation ───────────────────────────────────────────────
# MOFA needs: features (genes/peaks) as ROWS, samples as COLUMNS
# Check if your CSVs have samples as rows or columns

# Quick check: if ncol >> nrow, samples are probably columns (correct)
# If nrow >> ncol, samples are rows and we need to transpose
cat("\nRNA: ", nrow(rna), "rows x", ncol(rna), "cols\n")
cat("Neg: ", nrow(metabo_neg), "rows x", ncol(metabo_neg), "cols\n")
cat("Pos: ", nrow(metabo_pos), "rows x", ncol(metabo_pos), "cols\n")

# Transpose if needed (uncomment the relevant lines)
# If samples are ROWS and features are COLUMNS, transpose:
# rna <- t(rna)
# metabo_neg <- t(metabo_neg)
# metabo_pos <- t(metabo_pos)

# After this step: rows = features, columns = samples

# ── PART 5: Sample matching ─────────────────────────────────────────────────

# Find samples present in ALL three datasets
shared_samples <- Reduce(intersect, list(
  colnames(rna),
  colnames(metabo_neg),
  colnames(metabo_pos)
))
cat("\nShared samples across all 3 views:", length(shared_samples), "\n")
# Should be ~142

# Subset to shared samples
rna_shared <- as.matrix(rna[, shared_samples])
neg_shared <- as.matrix(metabo_neg[, shared_samples])
pos_shared <- as.matrix(metabo_pos[, shared_samples])

cat("RNA shared:", dim(rna_shared), "\n")
cat("Neg shared:", dim(neg_shared), "\n")
cat("Pos shared:", dim(pos_shared), "\n")

# ── PART 6: Create MOFA object ──────────────────────────────────────────────

# Put matrices into a named list
data_list <- list(
  RNA = rna_shared,
  Metabo_Neg = neg_shared,
  Metabo_Pos = pos_shared
)

# Create MOFA object
mofa_obj <- create_mofa(data_list)
print(mofa_obj)

# ── PART 7: Set MOFA options ───────────────────────────────────────────────

# Data options
data_opts <- get_default_data_options(mofa_obj)
data_opts$scale_views <- TRUE  # Scale each view to unit variance (important!)
cat("\nData options:\n")
print(data_opts)

# Model options
model_opts <- get_default_model_options(mofa_obj)
model_opts$num_factors <- 10
cat("\nModel options:\n")
print(model_opts)

# Training options
train_opts <- get_default_training_options(mofa_obj)
train_opts$seed <- 42
train_opts$convergence_mode <- "medium"
cat("\nTraining options:\n")
print(train_opts)

# Apply options
mofa_obj <- prepare_mofa(
  object = mofa_obj,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# ── PART 8: Run MOFA ───────────────────────────────────────────────────────

# This is the main step — should take 5-15 minutes on a laptop
cat("\n=== Running MOFA2 (this may take 5-15 minutes) ===\n")
mofa_trained <- run_mofa(mofa_obj, outfile = "mofa_model_corrected.hdf5",
                          use_basilisk = TRUE)  # use_basilisk handles Python env

cat("\n=== MOFA2 training complete! ===\n")

# ── PART 9: Variance decomposition ─────────────────────────────────────────

# This is the KEY diagnostic — how much variance each factor explains per view
var_explained <- plot_variance_explained(mofa_trained, plot_total = TRUE)
print(var_explained[[1]])  # Per-factor heatmap
print(var_explained[[2]])  # Total per view

# Save the heatmap
ggsave("mofa_variance_explained.png", var_explained[[1]],
       width = 8, height = 5, dpi = 300)
ggsave("mofa_variance_total.png", var_explained[[2]],
       width = 6, height = 4, dpi = 300)

# Get the actual numbers
r2 <- get_variance_explained(mofa_trained)
cat("\nVariance explained per factor per view:\n")
print(r2$r2_per_factor)

# ── PART 10: Load metadata ─────────────────────────────────────────────────

# Load sample metadata
meta <- read.csv("sample_sheet.csv", stringsAsFactors = FALSE)

# Match metadata to MOFA samples
mofa_meta <- data.frame(
  sample = samples_names(mofa_trained)[[1]],
  stringsAsFactors = FALSE
)

# Merge — adjust column names if your sample_sheet uses different headers
# Common columns: SampleID, Site, Description (which contains concentration)
mofa_meta <- merge(mofa_meta, meta, by.x = "sample", by.y = "SampleID", all.x = TRUE)

# Extract concentration from Description or create it manually
# Adjust this based on your actual metadata columns
# Example: if you have a "Description" column with "Control", "1x", "10x"
# mofa_meta$Concentration <- mofa_meta$Description

# If your metadata has a Site column:
# mofa_meta$Site <- mofa_meta$Site

# Add metadata to MOFA object
samples_metadata(mofa_trained) <- mofa_meta

# ── PART 11: Factor scatter plots ───────────────────────────────────────────

# Factor 1 vs Factor 2 coloured by concentration
p_conc <- plot_factors(mofa_trained, factors = c(1, 2), color_by = "Description") +
  ggtitle("MOFA2 — Factor 1 vs Factor 2 (by Concentration)") +
  theme_minimal(base_size = 14)
ggsave("mofa_scatter_concentration.png", p_conc, width = 8, height = 6, dpi = 300)

# Factor 1 vs Factor 2 coloured by site
p_site <- plot_factors(mofa_trained, factors = c(1, 2), color_by = "Site") +
  ggtitle("MOFA2 — Factor 1 vs Factor 2 (by Site)") +
  theme_minimal(base_size = 14)
ggsave("mofa_scatter_site.png", p_site, width = 8, height = 6, dpi = 300)

# ── PART 12: Factor 1 boxplot by site ──────────────────────────────────────

# Extract factor scores
factors_df <- as.data.frame(get_factors(mofa_trained)[[1]])
factors_df$sample <- rownames(factors_df)
factors_df <- merge(factors_df, mofa_meta, by = "sample")

# Boxplot of Factor 1 across sites
p_box <- ggplot(factors_df, aes(x = Site, y = Factor1, fill = Site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Factor 1 Scores across Danube Sites",
       x = "Site (D01=upstream → D12=downstream)",
       y = "Factor 1 score") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("mofa_factor1_by_site.png", p_box, width = 10, height = 6, dpi = 300)

# ── PART 13: Top feature loadings (THE MISSING PIECE) ──────────────────────

# Top RNA loadings for Factor 1
rna_weights_f1 <- plot_top_weights(mofa_trained, view = "RNA", factor = 1, nfeatures = 20)
ggsave("mofa_top_rna_loadings_F1.png", rna_weights_f1, width = 8, height = 6, dpi = 300)

# Top RNA loadings for Factor 3 (if it's a shared factor)
rna_weights_f3 <- plot_top_weights(mofa_trained, view = "RNA", factor = 3, nfeatures = 20)
ggsave("mofa_top_rna_loadings_F3.png", rna_weights_f3, width = 8, height = 6, dpi = 300)

# Top Negative metabolomics loadings for Factor 2
neg_weights_f2 <- plot_top_weights(mofa_trained, view = "Metabo_Neg", factor = 2, nfeatures = 20)
ggsave("mofa_top_neg_loadings_F2.png", neg_weights_f2, width = 8, height = 6, dpi = 300)

# Top Positive metabolomics loadings for Factor 1 (should now show signal!)
pos_weights_f1 <- plot_top_weights(mofa_trained, view = "Metabo_Pos", factor = 1, nfeatures = 20)
ggsave("mofa_top_pos_loadings_F1.png", pos_weights_f1, width = 8, height = 6, dpi = 300)

# ── PART 14: Extract weights as tables for annotation ──────────────────────

# Get all weights (loadings) as a data frame
weights <- get_weights(mofa_trained, as.data.frame = TRUE)

# Top 50 RNA genes driving Factor 1
top_rna_f1 <- weights %>%
  filter(view == "RNA", factor == "Factor1") %>%
  arrange(desc(abs(value))) %>%
  head(50)
write.csv(top_rna_f1, "top50_rna_genes_factor1.csv", row.names = FALSE)

# Top 50 negative metabolomics features driving Factor 2
top_neg_f2 <- weights %>%
  filter(view == "Metabo_Neg", factor == "Factor2") %>%
  arrange(desc(abs(value))) %>%
  head(50)
write.csv(top_neg_f2, "top50_neg_metabo_factor2.csv", row.names = FALSE)

# Top 50 positive metabolomics features driving Factor 1
top_pos_f1 <- weights %>%
  filter(view == "Metabo_Pos", factor == "Factor1") %>%
  arrange(desc(abs(value))) %>%
  head(50)
write.csv(top_pos_f1, "top50_pos_metabo_factor1.csv", row.names = FALSE)

cat("\nTop loading tables saved.\n")

# ── PART 15: Map gene IDs to human-readable names ──────────────────────────

# Load ortholog mapping (from Rayane)
orthologs <- read.csv("ortholog_mapping_anchored.csv", stringsAsFactors = FALSE)
cat("Ortholog mapping:", nrow(orthologs), "rows\n")

# Merge top RNA genes with ortholog names
# Adjust column names based on your actual ortholog file
# The 'feature' column from MOFA weights should match a Daphnia gene ID column
top_rna_f1_annotated <- merge(
  top_rna_f1,
  orthologs,
  by.x = "feature",
  by.y = names(orthologs)[1],  # adjust if the Dma gene ID column has a different name
  all.x = TRUE
)
write.csv(top_rna_f1_annotated, "top50_rna_genes_factor1_annotated.csv", row.names = FALSE)
cat("Annotated gene table saved.\n")

# ── PART 16: WGCNA cross-validation heatmap ────────────────────────────────

# Load WGCNA module eigengenes
eigengenes <- read.csv("module_eigengenes.csv", row.names = 1, check.names = FALSE)

# Get MOFA factor scores
factor_scores <- as.data.frame(get_factors(mofa_trained)[[1]])

# Match samples between MOFA and WGCNA
shared_wgcna <- intersect(rownames(factor_scores), rownames(eigengenes))
cat("\nShared samples (MOFA ∩ WGCNA):", length(shared_wgcna), "\n")

factor_scores_shared <- factor_scores[shared_wgcna, ]
eigengenes_shared <- eigengenes[shared_wgcna, ]

# Compute Pearson correlation matrix
cor_matrix <- cor(factor_scores_shared, eigengenes_shared, use = "pairwise.complete.obs")

# Plot heatmap
png("mofa_wgcna_integration.png", width = 10, height = 7, units = "in", res = 300)
pheatmap(
  cor_matrix,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 8,
  fontsize = 10,
  main = "MOFA2 Factors vs WGCNA Module Eigengenes\n(Pearson correlation)",
  angle_col = 45
)
dev.off()

cat("\nWGCNA integration heatmap saved.\n")

# ── PART 17: Summary statistics ─────────────────────────────────────────────

cat("\n============================================\n")
cat("MOFA2 ANALYSIS COMPLETE — SUMMARY\n")
cat("============================================\n")
cat("Samples used:", length(shared_samples), "\n")
cat("Views: RNA (", nrow(rna_shared), "features),",
    "Metabo_Neg (", nrow(neg_shared), "features),",
    "Metabo_Pos (", nrow(pos_shared), "features)\n")
cat("Factors: 10\n")
cat("\nVariance explained (R²) per factor per view:\n")
print(round(r2$r2_per_factor[[1]], 2))
cat("\nTotal variance explained per view:\n")
print(round(r2$r2_total[[1]], 2))
cat("\nFiles generated:\n")
cat("  mofa_model_corrected.hdf5 — trained model\n")
cat("  mofa_variance_explained.png — variance heatmap\n")
cat("  mofa_variance_total.png — total variance per view\n")
cat("  mofa_scatter_concentration.png — F1 vs F2 by concentration\n")
cat("  mofa_scatter_site.png — F1 vs F2 by site\n")
cat("  mofa_factor1_by_site.png — Factor 1 boxplot\n")
cat("  mofa_top_rna_loadings_F1.png — top RNA genes on Factor 1\n")
cat("  mofa_top_rna_loadings_F3.png — top RNA genes on Factor 3\n")
cat("  mofa_top_neg_loadings_F2.png — top neg metabo on Factor 2\n")
cat("  mofa_top_pos_loadings_F1.png — top pos metabo on Factor 1\n")
cat("  top50_rna_genes_factor1.csv — gene IDs\n")
cat("  top50_rna_genes_factor1_annotated.csv — with human orthologs\n")
cat("  top50_neg_metabo_factor2.csv — peak IDs\n")
cat("  top50_pos_metabo_factor1.csv — peak IDs\n")
cat("  mofa_wgcna_integration.png — cross-validation heatmap\n")
cat("============================================\n")
