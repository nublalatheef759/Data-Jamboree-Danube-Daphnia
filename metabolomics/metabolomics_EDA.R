knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.align = "center"
)

# ============================================================================
# Load Libraries and Data ===================================================
# ============================================================================

library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

# Set a clean colour palette for sites (12 Danube sites + control)
site_colors <- c(
  "CK"  = "#444444",   # controls = dark grey
  "D01" = "#1f77b4", "D02" = "#aec7e8", "D03" = "#ff7f0e",
  "D04" = "#ffbb78", "D05" = "#2ca02c", "D06" = "#98df8a",
  "D07" = "#d62728", "D08" = "#ff9896", "D09" = "#9467bd",
  "D10" = "#c5b0d5", "D11" = "#8c564b", "D12" = "#e377c2"
)

# ---- CHANGE THIS PATH to your actual file location ----
pos <- read.csv("polar_pos_pqn_imputed_glog.csv", row.names = 1, check.names = FALSE)

# The first column after row names is 'mz' (mass-to-charge ratio of each feature)
# All remaining columns are samples
mz_values  <- pos$mz
pos_mat    <- pos[ , colnames(pos) != "mz"]   # features x samples matrix

cat("Features (metabolite peaks):", nrow(pos_mat), "\n")
cat("Samples:                    ", ncol(pos_mat), "\n")
cat("Value range:                ", round(min(pos_mat), 3), "to", round(max(pos_mat), 3), "\n")

# ============================================================================
# Sample Metadata ===========================================================
# ============================================================================

# Parse sample names into structured metadata
samples <- colnames(pos_mat)

parse_sample <- function(s) {
  if (grepl("^CK", s)) {
    data.frame(sample = s, type = "Control", site = "CK",
               batch = substr(s, 3, 3), rep = substr(s, 4, 4),
               site_label = "CK (Control)", stringsAsFactors = FALSE)
  } else {
    data.frame(sample = s, type = "Danube", site = substr(s, 1, 3),
               batch = substr(s, 4, 4), rep = substr(s, 5, 5),
               site_label = substr(s, 1, 3), stringsAsFactors = FALSE)
  }
}

meta <- do.call(rbind, lapply(samples, parse_sample))
rownames(meta) <- meta$sample

cat("Sample breakdown:\n")
print(table(meta$site))

# ============================================================================
# Data Completeness Check ===================================================
# ============================================================================

# Missing values
n_missing <- sum(is.na(pos_mat))
cat("Missing values: ", n_missing, "\n")
cat("The data is already imputed — so this should be 0.\n")

# Per-sample total feature count (should be identical since fully imputed)
per_sample_nonmissing <- colSums(!is.na(pos_mat))
cat("\nPer-sample feature counts (should all be", nrow(pos_mat), "):\n")
print(summary(per_sample_nonmissing))

# Compare sample counts between pos and neg mode (run once with both loaded)
# Uncomment when running jointly:
# neg <- read.csv("polar_neg_pqn_imputed_glog.csv", row.names = 1, check.names = FALSE)
# neg_samples <- setdiff(colnames(neg), "mz")
# pos_only <- setdiff(colnames(pos_mat), neg_samples)
# neg_only <- setdiff(neg_samples, colnames(pos_mat))
# cat("Samples in POS but not NEG:", paste(pos_only, collapse=", "), "\n")
# cat("Samples in NEG but not POS:", paste(neg_only, collapse=", "), "\n")

# ============================================================================
# Per-Sample QC — Spotting Outliers =========================================
# ============================================================================

# Distribution of per-sample means ------------------------------------------

sample_means <- colMeans(pos_mat)
sample_sds   <- apply(pos_mat, 2, sd)

df_qc <- data.frame(
  sample = names(sample_means),
  mean   = sample_means,
  sd     = sample_sds,
  site   = meta[names(sample_means), "site"],
  batch  = meta[names(sample_means), "batch"]
)

grand_mean <- mean(df_qc$mean)
grand_sd   <- sd(df_qc$mean)

df_qc$outlier <- abs(df_qc$mean - grand_mean) > 2 * grand_sd

ggplot(df_qc, aes(x = reorder(sample, mean), y = mean,
                  fill = site, shape = batch)) +
  geom_col(aes(fill = site), alpha = 0.8) +
  geom_hline(yintercept = grand_mean + 2 * grand_sd,
             linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = grand_mean - 2 * grand_sd,
             linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_fill_manual(values = site_colors) +
  labs(title = "Per-sample mean intensity (positive mode)",
       subtitle = "Red dashed lines = ±2 SD from grand mean. Bars outside = potential outliers.",
       x = "Sample", y = "Mean glog intensity") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position = "right")

cat("Flagged outlier samples (mean > ±2 SD):\n")
print(df_qc[df_qc$outlier, c("sample", "mean", "sd", "site", "batch")])

# Boxplots per sample — are distributions comparable? -----------------------

# Take a random subset for readability (every 3rd sample)
idx <- seq(1, ncol(pos_mat), by = 3)
subset_mat <- pos_mat[ , idx]

df_box <- melt(as.data.frame(t(subset_mat)))
colnames(df_box) <- c("sample", "intensity")
df_box$site <- meta[as.character(df_box$sample), "site"]

ggplot(df_box, aes(x = sample, y = intensity, fill = site)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4) +
  scale_fill_manual(values = site_colors) +
  labs(title = "Intensity distribution per sample (every 3rd sample shown)",
       subtitle = "After PQN normalisation — boxes should be roughly aligned.",
       x = NULL, y = "glog intensity") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))

# ============================================================================
# Feature-Level QC ==========================================================
# ============================================================================

feat_var    <- apply(pos_mat, 1, var)
feat_mean   <- apply(pos_mat, 1, mean)
feat_cv     <- apply(pos_mat, 1, sd) / feat_mean  # coefficient of variation

cat("Feature variance summary:\n")
print(summary(feat_var))

# Near-zero variance features
nzv <- sum(feat_var < 0.01)
cat("\nNear-zero variance features (var < 0.01):", nzv, "\n")
cat("These features carry no information — consider removing for analysis.\n")

# Plot variance distribution
df_var <- data.frame(variance = feat_var, mean_intensity = feat_mean)

ggplot(df_var, aes(x = variance)) +
  geom_histogram(bins = 80, fill = "#2c7bb6", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0.01, color = "red", linetype = "dashed") +
  labs(title = "Distribution of feature variance across all samples",
       subtitle = "Features left of red line (var < 0.01) are near-constant — biologically uninformative.",
       x = "Variance", y = "Count") +
  theme_bw(base_size = 12)

ggplot(df_var, aes(x = mean_intensity, y = variance)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#2c7bb6") +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 0.8) +
  labs(title = "Mean–Variance relationship",
       subtitle = "In glog-transformed data this should be roughly flat. A strong trend = transformation didn't fully stabilise variance.",
       x = "Mean intensity", y = "Variance") +
  theme_bw(base_size = 12)

# ============================================================================
# Annotation Coverage =======================================================
# ============================================================================

pos_kegg <- read.table("polar_pos_pkl_to_kegg_annotations.tsv",
                        sep = "\t", header = FALSE,
                        col.names = c("feature_id", "mz", "kegg_ids"),
                        stringsAsFactors = FALSE)
pos_hmdb <- read.table("polar_pos_pkl_to_hmdb_annotations.tsv",
                        sep = "\t", header = FALSE,
                        col.names = c("feature_id", "mz", "hmdb_ids"),
                        stringsAsFactors = FALSE)

n_total <- nrow(pos_kegg)
n_kegg  <- sum(!is.na(pos_kegg$kegg_ids) & pos_kegg$kegg_ids != "")
n_hmdb  <- sum(!is.na(pos_hmdb$hmdb_ids) & pos_hmdb$hmdb_ids != "")
n_both  <- sum((!is.na(pos_kegg$kegg_ids) & pos_kegg$kegg_ids != "") &
               (!is.na(pos_hmdb$hmdb_ids) & pos_hmdb$hmdb_ids != ""))
n_none  <- sum((is.na(pos_kegg$kegg_ids) | pos_kegg$kegg_ids == "") &
               (is.na(pos_hmdb$hmdb_ids) | pos_hmdb$hmdb_ids == ""))

df_ann <- data.frame(
  category = c("KEGG only", "HMDB only", "Both KEGG+HMDB", "Unannotated"),
  count    = c(n_kegg - n_both, n_hmdb - n_both, n_both, n_none)
)

ggplot(df_ann, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Annotation coverage — positive mode (1,285 features)",
       fill = "Annotation type") +
  theme_void(base_size = 13)

cat("\nAnnotation breakdown:\n")
print(df_ann)
cat("\nNote: ~34% KEGG, ~32% HMDB annotation rate is typical for untargeted metabolomics.\n")
cat("Unannotated features are NOT useless — they still carry biological signal.\n")
cat("We just can't assign them a pathway name yet.\n")

# ============================================================================
# PCA — The Big Picture =====================================================
# ============================================================================

# PCA on samples (transpose: rows=samples, cols=features)
pca_res  <- prcomp(t(pos_mat), scale. = TRUE, center = TRUE)
pca_var  <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- as.data.frame(pca_res$x[ , 1:3])
pca_df$sample <- rownames(pca_df)
pca_df$site   <- meta[rownames(pca_df), "site"]
pca_df$batch  <- meta[rownames(pca_df), "batch"]
pca_df$type   <- meta[rownames(pca_df), "type"]

# PC1 vs PC2 coloured by site
ggplot(pca_df, aes(x = PC1, y = PC2, color = site, shape = batch)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(values = site_colors) +
  stat_ellipse(aes(group = site), level = 0.75, linewidth = 0.4, linetype = "dashed") +
  labs(title = "PCA — Positive mode metabolomics",
       subtitle = paste0("PC1: ", pca_var[1], "% variance | PC2: ", pca_var[2], "% variance | Coloured by site"),
       x = paste0("PC1 (", pca_var[1], "%)"),
       y = paste0("PC2 (", pca_var[2], "%)")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

# PC2 vs PC3 — sometimes reveals a different grouping
ggplot(pca_df, aes(x = PC2, y = PC3, color = site, shape = batch)) +
  geom_point(size = 3.5, alpha = 0.85) +
  scale_color_manual(values = site_colors) +
  labs(title = "PCA — PC2 vs PC3",
       subtitle = paste0("PC2: ", pca_var[2], "% | PC3: ", pca_var[3], "%"),
       x = paste0("PC2 (", pca_var[2], "%)"),
       y = paste0("PC3 (", pca_var[3], "%)")) +
  theme_bw(base_size = 12)

# Scree plot — how many PCs carry meaningful variance?
scree_df <- data.frame(PC = 1:20, variance = pca_var[1:20])

ggplot(scree_df, aes(x = PC, y = variance)) +
  geom_col(fill = "#2c7bb6", alpha = 0.8) +
  geom_line(color = "black", linewidth = 0.6) +
  geom_point(size = 2) +
  labs(title = "Scree plot — variance explained per PC",
       subtitle = "The 'elbow' tells you how many PCs contain real signal vs noise.",
       x = "Principal Component", y = "Variance explained (%)") +
  theme_bw(base_size = 12)

# PCA: Control vs Danube split ----------------------------------------------

ggplot(pca_df, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = c("Control" = "#333333", "Danube" = "#d73027")) +
  labs(title = "PCA — Control vs Danube samples",
       subtitle = "If controls cluster separately, pollution has a clear metabolic signature.",
       x = paste0("PC1 (", pca_var[1], "%)"),
       y = paste0("PC2 (", pca_var[2], "%)")) +
  theme_bw(base_size = 12)

# ============================================================================
# Outlier Deep-Dive =========================================================
# ============================================================================

# Mahalanobis distance from centroid in PC1-PC2 space
pc_scores <- pca_df[ , c("PC1", "PC2")]
centroid  <- colMeans(pc_scores)
cov_mat   <- cov(pc_scores)
mah_dist  <- mahalanobis(pc_scores, centroid, cov_mat)

pca_df$mahal  <- mah_dist
pca_df$pca_outlier <- mah_dist > qchisq(0.975, df = 2)  # chi-sq threshold

ggplot(pca_df, aes(x = PC1, y = PC2, color = pca_outlier, label = sample)) +
  geom_point(size = 3, alpha = 0.85) +
  ggrepel::geom_text_repel(data = pca_df[pca_df$pca_outlier, ],
                            size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("FALSE" = "#888888", "TRUE" = "#d73027"),
                     labels = c("Normal", "PCA outlier")) +
  labs(title = "PCA outlier detection (Mahalanobis distance)",
       subtitle = "Red points are >97.5th percentile distance from the PCA centroid.",
       color = "Status") +
  theme_bw(base_size = 12)

cat("PCA outliers (Mahalanobis > chi-sq 97.5%):\n")
outlier_tab <- pca_df[pca_df$pca_outlier, c("sample", "site", "batch", "PC1", "PC2", "mahal")]
outlier_tab <- outlier_tab[order(-outlier_tab$mahal), ]
print(round(outlier_tab, 3))

# ============================================================================
# Within-Site Reproducibility ===============================================
# ============================================================================

sites_danube <- unique(meta$site[meta$type == "Danube"])

within_site_cv <- sapply(sites_danube, function(s) {
  s_cols <- meta$sample[meta$site == s]
  s_cols <- s_cols[s_cols %in% colnames(pos_mat)]
  sub_mat <- pos_mat[ , s_cols, drop = FALSE]
  cv_vals <- apply(sub_mat, 1, function(x) sd(x) / mean(x))
  median(cv_vals, na.rm = TRUE)
})

df_cv <- data.frame(site = names(within_site_cv), median_cv = within_site_cv)

ggplot(df_cv, aes(x = site, y = median_cv, fill = site)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(values = site_colors) +
  labs(title = "Within-site variability (median CV across features)",
       subtitle = "Lower = replicates are more consistent. Spikes may reflect biological stress OR technical issues.",
       x = "Site", y = "Median coefficient of variation") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# ============================================================================
# Heatmap — Full Feature × Sample Matrix ====================================
# ============================================================================

# Use top 100 most variable features for clarity
top_var_idx <- order(feat_var, decreasing = TRUE)[1:100]
heat_mat    <- as.matrix(pos_mat[top_var_idx, ])

# Sample annotation bar
ann_col <- data.frame(
  Site  = meta[colnames(heat_mat), "site"],
  Batch = meta[colnames(heat_mat), "batch"],
  row.names = colnames(heat_mat)
)

ann_colors <- list(
  Site  = site_colors,
  Batch = c("A" = "#e41a1c", "B" = "#377eb8")
)

pheatmap(
  heat_mat,
  annotation_col   = ann_col,
  annotation_colors = ann_colors,
  show_rownames    = FALSE,
  show_colnames    = FALSE,
  scale            = "row",          # Z-score per feature — makes patterns visible
  clustering_method = "ward.D2",
  color            = colorRampPalette(c("#313695","#f7f7f7","#a50026"))(100),
  main             = "Top 100 most variable features × all samples\n(Z-scored by row, clustered by Ward's D2)"
)

# ============================================================================
# m/z Distribution ==========================================================
# ============================================================================

df_mz <- data.frame(
  mz          = mz_values,
  has_kegg    = !is.na(pos_kegg$kegg_ids) & pos_kegg$kegg_ids != "",
  has_hmdb    = !is.na(pos_hmdb$hmdb_ids) & pos_hmdb$hmdb_ids != ""
)

ggplot(df_mz, aes(x = mz, fill = has_kegg)) +
  geom_histogram(bins = 60, alpha = 0.8, position = "stack") +
  scale_fill_manual(values = c("FALSE" = "#aaaaaa", "TRUE" = "#2c7bb6"),
                    labels = c("Unannotated", "KEGG annotated")) +
  labs(title = "m/z distribution of features — positive mode",
       subtitle = "Annotated features shown in blue. Coverage drops at high m/z (larger molecules are harder to match).",
       x = "m/z (mass-to-charge ratio)", y = "Feature count", fill = "") +
  theme_bw(base_size = 12)

# ============================================================================
# Control vs Danube Comparison ==============================================
# ============================================================================

ctrl_cols   <- meta$sample[meta$type == "Control"]
danube_cols <- meta$sample[meta$type == "Danube"]

ctrl_means   <- rowMeans(pos_mat[ , ctrl_cols])
danube_means <- rowMeans(pos_mat[ , danube_cols])

fold_change  <- danube_means - ctrl_means   # log-scale difference = fold change

df_fc <- data.frame(
  feature     = rownames(pos_mat),
  mz          = mz_values,
  fold_change = fold_change,
  mean_ctrl   = ctrl_means
)

ggplot(df_fc, aes(x = mean_ctrl, y = fold_change)) +
  geom_point(alpha = 0.3, size = 0.7, color = "#555555") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(-0.5, 0.5), color = "red",
             linetype = "dashed", linewidth = 0.5) +
  labs(title = "MA plot — Danube vs Control (all sites combined)",
       subtitle = "Y-axis: mean Danube minus mean Control (in glog space = log fold change).\nRed dashed = ±0.5 log fold change. Points above = higher in Danube.",
       x = "Mean intensity (control)", y = "Log fold change (Danube − Control)") +
  theme_bw(base_size = 12)

# ============================================================================
# Site-Level Mean Profiles ==================================================
# ============================================================================

# Mean intensity per site for the top 50 variable features
top50_idx <- order(feat_var, decreasing = TRUE)[1:50]
top50_mat <- pos_mat[top50_idx, ]

site_means <- sapply(c("CK", paste0("D", sprintf("%02d", 1:12))), function(s) {
  s_cols <- meta$sample[meta$site == s]
  s_cols <- s_cols[s_cols %in% colnames(top50_mat)]
  if (length(s_cols) == 0) return(rep(NA, nrow(top50_mat)))
  rowMeans(top50_mat[ , s_cols, drop = FALSE])
})

df_site <- melt(site_means)
colnames(df_site) <- c("feature_idx", "site", "mean_intensity")

ggplot(df_site, aes(x = site, y = mean_intensity, fill = site)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
  scale_fill_manual(values = c(site_colors,
                                "D01"="#1f77b4","D02"="#aec7e8","D03"="#ff7f0e",
                                "D04"="#ffbb78","D05"="#2ca02c","D06"="#98df8a",
                                "D07"="#d62728","D08"="#ff9896","D09"="#9467bd",
                                "D10"="#c5b0d5","D11"="#8c564b","D12"="#e377c2")) +
  labs(title = "Distribution of top-50 feature intensities by site",
       subtitle = "Each box = one site's mean intensity across the 50 most variable features.\nSite ordering = approximate Danube river position (D01 upstream → D12 downstream).",
       x = "Site", y = "Mean intensity (top 50 features)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# ============================================================================
# Annotation Summary Table ==================================================
# ============================================================================

# Build a reference table of your annotated features
annot_tab <- merge(pos_kegg, pos_hmdb, by = "feature_id", suffixes = c("_kegg", "_hmdb"))
annot_tab <- annot_tab[ , c("feature_id", "mz_kegg", "kegg_ids", "hmdb_ids")]
colnames(annot_tab) <- c("feature_id", "mz", "kegg_ids", "hmdb_ids")

has_any <- !is.na(annot_tab$kegg_ids) | !is.na(annot_tab$hmdb_ids)
annot_tab <- annot_tab[has_any, ]

cat("Features with at least one annotation:", nrow(annot_tab), "\n")
cat("Sample of annotated features:\n")
print(head(annot_tab[!is.na(annot_tab$kegg_ids), ], 10))

# ============================================================================
# EDA Summary — What We Found ===============================================
# ============================================================================

summary_df <- data.frame(
  Item = c(
    "Total features (positive mode)",
    "Total samples",
    "Missing values",
    "Near-zero variance features",
    "KEGG annotated features",
    "HMDB annotated features",
    "m/z range",
    "Value range (glog transformed)",
    "Samples in POS but not NEG",
    "Samples in NEG but not POS"
  ),
  Value = c(
    "1,285",
    "149",
    "0 (fully imputed)",
    "0 (positive mode)",
    "437 (34.0%)",
    "409 (31.8%)",
    "52.5 – 617.2",
    "7.98 – 15.39",
    "D06A2, D06B3, D07A3, D08B1, D10B4",
    "D01B4"
  )
)
knitr::kable(summary_df, caption = "EDA Summary — Positive Mode Metabolomics")

# What to tell your teammates -----------------------------------------------

# ============================================================================
# Session Info ==============================================================
# ============================================================================

sessionInfo()
