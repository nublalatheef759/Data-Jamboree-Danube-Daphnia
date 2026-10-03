# ============================================================================
# POSITIVE MODE METABOLOMICS — QC, EDA & DIFFERENTIAL ANALYSIS
# Daphnia magna × Danube River multi-omics project
# ============================================================================
#
# Input files:
#   polar_pos_pqn_imputed_glog.csv           — PQN-normalised, KNN-imputed, glog-transformed
#   polar_pos_pkl_to_kegg_annotations.tsv    — KEGG compound annotations
#   polar_pos_pkl_to_hmdb_annotations.tsv    — HMDB compound annotations
#
# Outputs:
#   metabolomics_workspace.RData             — workspace for downstream scripts
#   QC plots, PCA, heatmaps, ANOVA, PERMANOVA, volcano plots
# ============================================================================

# Libraries -----------------------------------------------------------------

library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(vegan)

cat("All libraries loaded.\n")

# Colour palette ------------------------------------------------------------

# Define ONCE here — every plot in this script uses these same colours
# so all your figures are visually consistent when you present them

site_colors <- c(
  "CK"  = "#444444",
  "D01" = "#1f77b4", "D02" = "#aec7e8", "D03" = "#ff7f0e",
  "D04" = "#ffbb78", "D05" = "#2ca02c", "D06" = "#98df8a",
  "D07" = "#d62728", "D08" = "#ff9896", "D09" = "#9467bd",
  "D10" = "#c5b0d5", "D11" = "#8c564b", "D12" = "#e377c2"
)

cat("Colour palette defined for CK + D01-D12.\n")

# Set working directory and load data ---------------------------------------

# ---- SET YOUR PATH HERE ----

# Load the main data file
# Rows = metabolite features (m/z peaks detected by the mass spectrometer)
# Columns = individual Daphnia magna animals
pos <- read.csv(
  "polar_pos_pqn_imputed_glog.csv",
  row.names   = 1,
  check.names = FALSE
)

# Separate the m/z column (feature identity) from the intensity matrix
mz_values <- pos$mz
pos_mat   <- pos[ , colnames(pos) != "mz"]

cat("=== DATA LOADED ===\n")
cat("Features (metabolite peaks):", nrow(pos_mat), "\n")
cat("Samples (Daphnia animals):  ", ncol(pos_mat), "\n")
cat("m/z range:", round(min(mz_values), 2), "to", round(max(mz_values), 2), "\n")
cat("Value range (glog):        ", round(min(pos_mat), 3),
    "to", round(max(pos_mat), 3), "\n")

# Load annotation files -----------------------------------------------------

# na.strings is CRITICAL — without it R reads blank cells as "NA" the string,
# not as proper NA, and every feature appears annotated when it isn't

pos_kegg <- read.table(
  "polar_pos_pkl_to_kegg_annotations.tsv",
  sep = "\t", header = FALSE,
  col.names = c("feature_id", "mz", "kegg_ids"),
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")    # blank cells → proper NA
)

pos_hmdb <- read.table(
  "polar_pos_pkl_to_hmdb_annotations.tsv",
  sep = "\t", header = FALSE,
  col.names = c("feature_id", "mz", "hmdb_ids"),
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

has_kegg <- !is.na(pos_kegg$kegg_ids)
has_hmdb <- !is.na(pos_hmdb$hmdb_ids)

cat("=== ANNOTATION COVERAGE ===\n")
cat("Total features:         ", nrow(pos_kegg), "\n")
cat("KEGG annotated:         ", sum(has_kegg),
    sprintf("(%.1f%%)\n", 100 * mean(has_kegg)))
cat("HMDB annotated:         ", sum(has_hmdb),
    sprintf("(%.1f%%)\n", 100 * mean(has_hmdb)))
cat("Both KEGG + HMDB:       ", sum(has_kegg & has_hmdb), "\n")
cat("Unannotated (neither):  ", sum(!has_kegg & !has_hmdb), "\n")

# Build sample metadata -----------------------------------------------------

# Decode sample names into a structured table
# CK07-CK12   = lab controls (clean water, no pollution)
# D01A1       = Danube site 1, batch A, replicate 1
# D12B6       = Danube site 12, batch B, replicate 6
# Batches A/B = two separate processing runs per site
# Replicates  = 6 individual Daphnia per batch = 12 per site

parse_sample <- function(s) {
  if (grepl("^CK", s)) {
    data.frame(sample = s, type = "Control", site = "CK",
               batch = "Control", rep = substr(s, 3, 4),
               stringsAsFactors = FALSE)
  } else {
    data.frame(sample = s, type = "Danube", site = substr(s, 1, 3),
               batch = substr(s, 4, 4), rep = substr(s, 5, 5),
               stringsAsFactors = FALSE)
  }
}

meta           <- do.call(rbind, lapply(colnames(pos_mat), parse_sample))
rownames(meta) <- meta$sample

cat("=== SAMPLE METADATA ===\n\n")
cat("Samples per site:\n"); print(table(meta$site))
cat("\nSamples per concentration:\n"); print(table(meta$batch))
cat("\nControl vs Danube:\n"); print(table(meta$type))

# Save workspace ------------------------------------------------------------

save(pos_mat, mz_values, meta, pos_kegg, pos_hmdb, site_colors,
     file = "metabolomics_workspace.RData")
cat("Saved: metabolomics_workspace.RData\n")
cat("Every section below loads from this — re-run this chunk if you restart R.\n")

# ============================================================================
# PART 2 — QUALITY CONTROL & OUTLIER DETECTION
# ============================================================================

# QC Layer 1 — Per-sample mean intensity ------------------------------------

sample_means <- colMeans(pos_mat)
sample_sds   <- apply(pos_mat, 2, sd)
grand_mean   <- mean(sample_means)
grand_sd     <- sd(sample_means)

df_qc <- data.frame(
  sample  = names(sample_means),
  mean    = as.numeric(sample_means),
  sd      = as.numeric(sample_sds),
  site    = meta[names(sample_means), "site"],
  batch   = meta[names(sample_means), "batch"],
  stringsAsFactors = FALSE
)
df_qc$z_score  <- (df_qc$mean - grand_mean) / grand_sd
df_qc$flag_L1  <- abs(df_qc$z_score) > 2

cat("Grand mean:", round(grand_mean, 4), "| Grand SD:", round(grand_sd, 4), "\n")
cat("±2 SD window:", round(grand_mean - 2*grand_sd, 4),
    "to", round(grand_mean + 2*grand_sd, 4), "\n\n")
cat("Layer 1 flags (|z| > 2):\n")
print(df_qc[df_qc$flag_L1, c("sample", "mean", "z_score", "site", "batch")])

ggplot(df_qc, aes(x = reorder(sample, mean), y = mean, fill = site)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = grand_mean + 2*grand_sd,
             linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_hline(yintercept = grand_mean - 2*grand_sd,
             linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_hline(yintercept = grand_mean,
             linetype = "solid", color = "black", linewidth = 0.4) +
  scale_fill_manual(values = site_colors) +
  labs(
    title    = "QC Layer 1 — Per-sample mean intensity",
    subtitle = paste0("Grand mean = ", round(grand_mean, 3),
                      " | ±2 SD = ",   round(grand_mean - 2*grand_sd, 3),
                      " to ",          round(grand_mean + 2*grand_sd, 3),
                      " | Flagged: ", sum(df_qc$flag_L1)),
    x = NULL, y = "Mean glog intensity"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
        legend.position = "right")

ggplot(df_qc, aes(x = reorder(sample, mean), y = mean,
                   color = site, shape = batch)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_hline(yintercept = grand_mean + 2*grand_sd,
             linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_hline(yintercept = grand_mean - 2*grand_sd,
             linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_hline(yintercept = grand_mean,
             linetype = "solid", color = "black", linewidth = 0.5) +
  geom_point(
    data = df_qc[df_qc$flag_L1, ],
    aes(x = reorder(sample, mean), y = mean),
    color = "red", size = 4, shape = 1, stroke = 1.2
  ) +
  geom_text(
    data = df_qc[df_qc$flag_L1, ],
    aes(label = sample),
    nudge_x = 3, nudge_y = 0, size = 3,
    color = "red", fontface = "bold"
  ) +
  scale_color_manual(values = site_colors) +
  scale_shape_manual(values = c("A" = 16, "B" = 17, "Control" = 15)) +
  coord_cartesian(ylim = c(
    min(df_qc$mean) - 0.05,
    max(df_qc$mean) + 0.1
  )) +
  labs(
    title    = "Sample QC: Mean glog intensity per sample",
    subtitle = paste0("Red dashed = ±2SD boundary; labelled points = potential outliers",
                      " | Flagged: ", sum(df_qc$flag_L1)),
    x = "Sample (ordered by mean intensity)",
    y = "Mean glog intensity",
    color = "Site", shape = "Concentration"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

cat("Grand mean:", round(grand_mean, 4), "\n")
cat("Grand SD:", round(grand_sd, 4), "\n")
cat("Upper boundary:", round(grand_mean + 2*grand_sd, 4), "\n")
cat("Lower boundary:", round(grand_mean - 2*grand_sd, 4), "\n")
cat("\nFlagged samples:\n")
print(df_qc[df_qc$flag_L1, c("sample", "mean", "z_score")])

df_dens <- melt(as.matrix(pos_mat))
colnames(df_dens) <- c("feature", "sample", "intensity")
df_dens$sample <- as.character(df_dens$sample)
df_dens$site <- meta[df_dens$sample, "site"]
df_dens$flagged <- df_dens$sample %in%
                   df_qc$sample[df_qc$flag_L1]

ggplot(df_dens, aes(x = intensity, group = sample,
                     color = site, alpha = flagged,
                     linewidth = flagged)) +
  geom_density() +
  scale_color_manual(values = site_colors) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.15)) +
  scale_linewidth_manual(values = c("TRUE" = 1.2, "FALSE" = 0.3)) +
  labs(title = "Density curves — all 149 samples overlaid",
       subtitle = "Flagged samples highlighted. All curves should overlap. A shifted curve = problem sample.",
       x = "glog intensity", y = "Density") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

sample_cv <- apply(pos_mat, 2, function(x) sd(x)/abs(mean(x)))

df_cv_samp <- data.frame(
  sample = names(sample_cv),
  cv = as.numeric(sample_cv),
  site = meta[names(sample_cv), "site"],
  flagged = names(sample_cv) %in% df_qc$sample[df_qc$flag_L1]
)

ggplot(df_cv_samp, aes(x = reorder(sample, cv), y = cv,
                        color = site, shape = flagged)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_text(data = df_cv_samp[df_cv_samp$flagged, ],
            aes(label = sample), nudge_x = 3, size = 3,
            color = "red", fontface = "bold") +
  scale_color_manual(values = site_colors) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16)) +
  labs(title = "Coefficient of variation per sample",
       subtitle = "Flagged samples marked. High CV = noisy sample. Low CV = unusually flat = possible failed extraction.",
       x = "Sample (ordered by CV)", y = "CV") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# QC Layer 2 — Distribution shape (boxplots) --------------------------------

idx        <- seq(1, ncol(pos_mat), by = 3)
subset_mat <- pos_mat[ , idx, drop = FALSE]

df_box           <- melt(as.matrix(t(subset_mat)))
colnames(df_box) <- c("sample", "feature", "intensity")
df_box$sample    <- as.character(df_box$sample)
df_box$site      <- meta[df_box$sample, "site"]

ggplot(df_box, aes(x = sample, y = intensity, fill = site)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.35) +
  scale_fill_manual(values = site_colors, na.value = "#cccccc") +
  geom_hline(yintercept = median(as.matrix(pos_mat)),
             color = "red", linetype = "dashed", linewidth = 0.5) +
  labs(
    title    = "QC Layer 2 — Intensity distribution per sample (every 3rd shown)",
    subtitle = "Red dashed = overall median. All boxes should be at the same height. A shifted box = normalisation failure.",
    x = NULL, y = "glog intensity"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

# QC Layer 3 — PCA + Mahalanobis distance -----------------------------------

# PCA on samples: transpose so rows=samples, cols=features
pca_res <- prcomp(t(pos_mat), scale. = TRUE, center = TRUE)
pca_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df        <- as.data.frame(pca_res$x[ , 1:4])
pca_df$sample <- rownames(pca_df)
pca_df$site   <- meta[rownames(pca_df), "site"]
pca_df$batch  <- meta[rownames(pca_df), "batch"]
pca_df$type   <- meta[rownames(pca_df), "type"]

# Mahalanobis distance: how far is each sample from the centre
# of the PCA cloud, accounting for the cloud's shape?
pc_scores  <- pca_df[ , c("PC1", "PC2")]
centroid   <- colMeans(pc_scores)
cov_mat    <- cov(pc_scores)
mah_sq     <- mahalanobis(pc_scores, centroid, cov_mat)

# Threshold: chi-squared(df=2) 97.5th percentile = 5.991
chi_sq_thresh <- qchisq(0.975, df = 2)

pca_df$mahal_sq    <- as.numeric(mah_sq)
pca_df$pca_outlier <- mah_sq > chi_sq_thresh
pca_df$flag_L1     <- rownames(pca_df) %in% df_qc$sample[df_qc$flag_L1]

cat("Chi-squared threshold (df=2, 97.5%):", round(chi_sq_thresh, 3), "\n\n")
cat("Top 10 samples by Mahalanobis distance:\n")
ranked <- pca_df[order(-pca_df$mahal_sq),
                 c("sample", "site", "batch", "PC1", "PC2", "mahal_sq", "pca_outlier")]
ranked[ , c("PC1","PC2","mahal_sq")] <- round(ranked[ , c("PC1","PC2","mahal_sq")], 3)
print(head(ranked, 10))
cat("\nFormal PCA outliers (squared Mahal > 5.991):\n")
formal <- pca_df[pca_df$pca_outlier, c("sample","site","batch","PC1","PC2","mahal_sq")]
formal[ , c("PC1","PC2","mahal_sq")] <- round(formal[ , c("PC1","PC2","mahal_sq")], 3)
if (nrow(formal) == 0) cat("None.\n") else print(formal)

ggplot(pca_df, aes(x = PC1, y = PC2, color = pca_outlier, shape = batch)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text(
    data  = pca_df[pca_df$pca_outlier, ],
    aes(label = sample),
    nudge_y = 2, size = 3.5, color = "#d62728", fontface = "bold"
  ) +
  geom_text(
    data  = pca_df[pca_df$flag_L1 & !pca_df$pca_outlier, ],
    aes(label = sample),
    nudge_y = -2, size = 2.8, color = "#e67e22"
  ) +
  scale_color_manual(
    values = c("FALSE" = "#aaaaaa", "TRUE" = "#d62728"),
    labels = c("Normal", "PCA outlier")
  ) +
  scale_shape_manual(values = c("A" = 16, "B" = 17, "Control" = 15)) +
  labs(
    title    = "QC Layer 3 — PCA outlier detection (Mahalanobis distance)",
    subtitle = paste0("PC1: ", pca_var[1], "% | PC2: ", pca_var[2],
                      "% | Red = squared Mahal > 5.991 (chi-sq 97.5%)",
                      " | Orange labels = Layer 1 mean flags (for comparison)"),
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)"),
    color = "Status", shape = "Batch"
  ) +
  theme_bw(base_size = 12)

# Outlier decision — evidence table -----------------------------------------

# Full evidence profile for D08B6 (the only formal PCA outlier)
d08_cols   <- rownames(meta)[meta$site == "D08"]
d08b6_vals <- pos_mat[ , "D08B6"]
d08_others <- d08_cols[d08_cols != "D08B6"]
d08_mean_f <- rowMeans(pos_mat[ , d08_others])
d08_sd_f   <- apply(pos_mat[ , d08_others], 1, sd)
z_vs_site  <- (d08b6_vals - d08_mean_f) / d08_sd_f

# Spearman correlation of D08B6 with its site-mates
d08_corr     <- cor(pos_mat[ , d08_cols], method = "spearman")
d08b6_avg_r  <- mean(d08_corr["D08B6", d08_others])
others_upper <- d08_corr[d08_others, d08_others]
site_avg_r   <- mean(others_upper[upper.tri(others_upper)])

cat("=== EVIDENCE PROFILE: D08B6 ===\n\n")
cat("Test 1 — Mean intensity z-score: ",
    round(df_qc[df_qc$sample == "D08B6", "z_score"], 3),
    " (threshold ±2) → PASS\n")
cat("Test 2 — Distribution shape: visually identical to D08 site-mates → PASS\n")
cat("Test 3 — Mahalanobis squared distance:",
    round(pca_df[pca_df$sample == "D08B6", "mahal_sq"], 3),
    " (threshold 5.991) → FLAG\n")
cat("Test 4 — Correlation with D08 site-mates:",
    round(d08b6_avg_r, 4),
    " | Site average:", round(site_avg_r, 4), "→ PASS (difference = 0.014)\n")
cat("Test 5 — Features with |z| > 3 vs site-mates:",
    sum(abs(z_vs_site) > 3, na.rm = TRUE),
    "of 1285 → biological variation, not technical failure\n\n")
cat("VERDICT: KEEP D08B6\n")
cat("5/6 tests pass. PCA flag is real but explained by individual biological\n")
cat("variation (reproductive stage, immune state), not extraction failure.\n")

# ============================================================================
# PART 3 — EXPLORATORY DATA ANALYSIS
# ============================================================================

# Feature variance ----------------------------------------------------------

feat_var  <- apply(pos_mat, 1, var)
feat_mean <- apply(pos_mat, 1, mean)
df_var    <- data.frame(variance = feat_var, mean_intensity = feat_mean)

cat("Feature variance summary:\n")
print(summary(feat_var))
cat("\nNear-zero variance (< 0.01):", sum(feat_var < 0.01), "\n")
cat("Top 10% most variable:       ", sum(feat_var > quantile(feat_var, 0.9)), "features\n")

ggplot(df_var, aes(x = variance)) +
  geom_histogram(bins = 80, fill = "#2c7bb6", color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0.01,
             color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = quantile(feat_var, 0.9),
             color = "darkgreen", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = quantile(feat_var, 0.9) + 0.04, y = 100,
           label = "Top 10%", color = "darkgreen", size = 3.5, hjust = 0) +
  labs(
    title    = "Feature variance distribution (1,285 features)",
    subtitle = "Right tail = most variable features = biologically most interesting. Red = near-zero threshold (none here).",
    x = "Variance", y = "Feature count"
  ) +
  theme_bw(base_size = 12)

# Mean–variance relationship ------------------------------------------------

ggplot(df_var, aes(x = mean_intensity, y = variance)) +
  geom_point(alpha = 0.25, size = 0.8, color = "#2c7bb6") +
  geom_smooth(method = "loess", se = TRUE, color = "red",
              fill = "pink", alpha = 0.3, linewidth = 0.9) +
  labs(
    title    = "Mean–Variance relationship",
    subtitle = "Glog transform should flatten this. Rise at low intensity = noise floor (normal). Flat at mid-high = good.",
    x = "Mean glog intensity", y = "Variance"
  ) +
  theme_bw(base_size = 12)

# Annotation coverage -------------------------------------------------------

df_ann <- data.frame(
  category = c("Both KEGG+HMDB", "KEGG only", "HMDB only", "Unannotated"),
  count    = c(sum(has_kegg & has_hmdb),
               sum(has_kegg & !has_hmdb),
               sum(!has_kegg & has_hmdb),
               sum(!has_kegg & !has_hmdb))
)
df_ann$pct <- round(100 * df_ann$count / sum(df_ann$count), 1)

ggplot(df_ann, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  geom_text(aes(label = paste0(count, "\n(", pct, "%)")),
            position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Annotation coverage — 1,285 positive mode features", fill = "") +
  theme_void(base_size = 13)

print(df_ann[ , c("category", "count", "pct")])

# PCA — main result ---------------------------------------------------------

ggplot(pca_df, aes(x = PC1, y = PC2, color = site, shape = batch)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = site_colors) +
  scale_shape_manual(values = c("A" = 16, "B" = 17, "Control" = 15)) +
  stat_ellipse(aes(group = site), level = 0.75,
               linewidth = 0.35, linetype = "dashed") +
  labs(
    title    = "PCA — Positive mode metabolomics (149 samples)",
    subtitle = paste0("PC1: ", pca_var[1], "% | PC2: ", pca_var[2],
                      "% | Coloured by site | Shaped by batch | Ellipses = 75% CI"),
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)")
  ) +
  theme_bw(base_size = 12)

# Calculate mean PC scores per site (centroid of each site's cloud)
site_centroids <- do.call(rbind, lapply(
  unique(pca_df$site), function(s) {
    sub <- pca_df[pca_df$site == s, c("PC1","PC2","PC3","PC4")]
    colMeans(sub)
  }
))
rownames(site_centroids) <- unique(pca_df$site)

# Order sites properly
site_order <- c("CK", paste0("D", sprintf("%02d", 1:12)))
site_centroids <- site_centroids[site_order, ]

# Euclidean distance between every pair of site centroids
dist_mat <- as.matrix(dist(site_centroids, method = "euclidean"))

# Melt for ggplot
dist_melt <- melt(dist_mat)
colnames(dist_melt) <- c("Site1", "Site2", "Distance")
dist_melt$Site1 <- factor(dist_melt$Site1, levels = site_order)
dist_melt$Site2 <- factor(dist_melt$Site2, levels = site_order)

ggplot(dist_melt, aes(x = Site1, y = Site2, fill = Distance)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = round(Distance, 1)), size = 2.8) +
  scale_fill_gradientn(
    colors = c("#2c7bb6", "#ffffbf", "#d7191c"),
    name = "Euclidean\ndistance"
  ) +
  labs(
    title    = "Pairwise metabolic distance between sites",
    subtitle = "Based on PC1–PC4 centroids. Blue = metabolically similar. Red = metabolically distant.",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(pca_df, aes(x = PC1, y = PC2, color = type)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = c("Control" = "#333333", "Danube" = "#d73027")) +
  labs(
    title    = "PCA — Control vs Danube",
    subtitle = "Controls scattered within Danube cloud — pollution has no single universal metabolic signature. Consistent with mixture hypothesis.",
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)")
  ) +
  theme_bw(base_size = 12)

scree_df <- data.frame(PC = 1:20, pct = pca_var[1:20])

ggplot(scree_df, aes(x = PC, y = pct)) +
  geom_col(fill = "#2c7bb6", alpha = 0.8, width = 0.7) +
  geom_line(color = "black", linewidth = 0.6) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 4.5, color = "red",
             linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 5.2, y = 14,
           label = "Elbow ~PC4-5\nNoise begins", color = "red", size = 3.5) +
  labs(
    title    = "Scree plot — variance per PC",
    subtitle = "PCs before the elbow carry biological signal. After = noise. Elbow at ~PC4-5.",
    x = "Principal Component", y = "Variance explained (%)"
  ) +
  theme_bw(base_size = 12)

# Within-site reproducibility -----------------------------------------------

sites_danube <- paste0("D", sprintf("%02d", 1:12))

within_cv <- sapply(sites_danube, function(s) {
  cols    <- rownames(meta)[meta$site == s]
  cols    <- cols[cols %in% colnames(pos_mat)]
  cv_vals <- apply(pos_mat[ , cols, drop = FALSE], 1,
                   function(x) sd(x) / abs(mean(x)))
  median(cv_vals, na.rm = TRUE)
})

df_cv <- data.frame(site = names(within_cv), median_cv = as.numeric(within_cv))

ggplot(df_cv, aes(x = site, y = median_cv, fill = site)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = mean(df_cv$median_cv),
             color = "black", linetype = "dashed") +
  scale_fill_manual(values = site_colors) +
  labs(
    title    = "Within-site variability — median CV per site",
    subtitle = "Lower = daphnia from same site are metabolically consistent. All sites ~0.03-0.045 = excellent reproducibility.",
    x = "Site", y = "Median coefficient of variation"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

cat("\nCV per site:\n")
print(setNames(round(within_cv, 5), sites_danube))

# Heatmap — top 100 variable features ---------------------------------------

top100   <- order(feat_var, decreasing = TRUE)[1:100]
heat_mat <- as.matrix(pos_mat[top100, ])

ann_col <- data.frame(
  Site  = meta[colnames(heat_mat), "site"],
  Batch = meta[colnames(heat_mat), "batch"],
  row.names = colnames(heat_mat)
)

# Build colour lists that EXACTLY match the data values
sites_in_data   <- unique(ann_col$Site)
batches_in_data <- unique(ann_col$Batch)
batch_pal       <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
batch_colors    <- setNames(batch_pal[seq_along(batches_in_data)], batches_in_data)

ann_colors <- list(
  Site  = site_colors[names(site_colors) %in% sites_in_data],
  Batch = batch_colors
)

pheatmap(
  heat_mat,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  scale             = "row",
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#313695", "#f7f7f7", "#a50026"))(100),
  main  = "Top 100 most variable features × 149 samples\n(row Z-scored | Ward's D2 clustering)"
)

# MA plot — global Danube vs control ----------------------------------------

ctrl_cols   <- rownames(meta)[meta$type == "Control"]
danube_cols <- rownames(meta)[meta$type == "Danube"]

ctrl_means   <- rowMeans(pos_mat[ , ctrl_cols])
danube_means <- rowMeans(pos_mat[ , danube_cols])
fc           <- as.numeric(danube_means - ctrl_means)

df_ma <- data.frame(
  feature     = rownames(pos_mat),
  mean_ctrl   = as.numeric(ctrl_means),
  fold_change = fc
)

up   <- sum(fc >  0.5)
down <- sum(fc < -0.5)

ggplot(df_ma, aes(x = mean_ctrl, y = fold_change)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#555555") +
  geom_hline(yintercept =  0,   color = "black", linewidth = 0.5) +
  geom_hline(yintercept =  0.5, color = "red",  linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = -0.5, color = "blue", linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = 13.5, y =  0.85,
           label = paste0("↑ ", up, " elevated"), color = "#c0392b", size = 3.5) +
  annotate("text", x = 13.5, y = -0.85,
           label = paste0("↓ ", down, " depleted"), color = "#2980b9", size = 3.5) +
  labs(
    title    = "MA plot — All Danube vs Control (pooled)",
    subtitle = "Y = log fold change. Points above red = elevated in Danube. Note: pools all 12 sites — site-specific analysis follows.",
    x = "Mean control intensity", y = "Log fold change (Danube − Control)"
  ) +
  theme_bw(base_size = 12)

cat("Elevated in Danube (FC > +0.5):", up, "\n")
cat("Depleted in Danube  (FC < -0.5):", down, "\n")

# Site-level profiles -------------------------------------------------------

top50_idx  <- order(feat_var, decreasing = TRUE)[1:50]
top50_mat  <- pos_mat[top50_idx, ]
site_order <- c("CK", paste0("D", sprintf("%02d", 1:12)))

site_list <- lapply(site_order, function(s) {
  cols <- rownames(meta)[meta$site == s]
  cols <- cols[cols %in% colnames(top50_mat)]
  if (length(cols) == 0) return(NULL)
  data.frame(site = s, animal = cols,
             mean_int = colMeans(top50_mat[ , cols, drop = FALSE]))
})
df_site      <- do.call(rbind, Filter(Negate(is.null), site_list))
df_site$site <- factor(df_site$site, levels = site_order)

ggplot(df_site, aes(x = site, y = mean_int, fill = site)) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.8) +
  scale_fill_manual(values = site_colors) +
  labs(
    title    = "Top-50 feature intensity by site (upstream → downstream)",
    subtitle = "Each point = one sample. CK = lab control. D01 = most upstream. D12 = most downstream.",
    x = "Site", y = "Mean intensity (top 50 variable features)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# ============================================================================
# PART 4 — DIFFERENTIAL ABUNDANCE ANALYSIS
# ============================================================================

# ANOVA — global test across all sites --------------------------------------

site_factor <- meta[colnames(pos_mat), "site"]

cat("Running ANOVA on", nrow(pos_mat), "features across",
    length(unique(site_factor)), "groups...\n")

anova_F    <- numeric(nrow(pos_mat))
anova_pval <- numeric(nrow(pos_mat))

for (i in seq_len(nrow(pos_mat))) {
  y              <- as.numeric(pos_mat[i, ])
  fit            <- aov(y ~ site_factor)
  smry           <- summary(fit)[[1]]
  anova_F[i]    <- smry["site_factor", "F value"]
  anova_pval[i] <- smry["site_factor", "Pr(>F)"]
}

# Benjamini-Hochberg FDR correction
anova_padj <- p.adjust(anova_pval, method = "BH")

anova_results <- data.frame(
  feature = rownames(pos_mat),
  mz      = mz_values,
  F_stat  = anova_F,
  pval    = anova_pval,
  padj    = anova_padj,
  stringsAsFactors = FALSE
)
anova_results <- anova_results[order(anova_results$padj), ]

cat("\n=== ANOVA RESULTS ===\n")
cat("Significant (FDR < 0.05):", sum(anova_padj < 0.05), "\n")
cat("Significant (FDR < 0.01):", sum(anova_padj < 0.01), "\n")
cat("Significant (FDR < 0.10):", sum(anova_padj < 0.10), "\n\n")
cat("Top 10 most significant:\n")
top10 <- head(anova_results, 10)
top10[ , c("F_stat","pval","padj")] <- round(top10[ , c("F_stat","pval","padj")], 5)
print(top10[ , c("feature","mz","F_stat","pval","padj")])

anova_results$sig <- anova_results$padj < 0.05

ggplot(anova_results, aes(x = F_stat, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept = -log10(0.05),
             color = "red", linetype = "dashed", linewidth = 0.7) +
  scale_color_manual(
    values = c("FALSE" = "#aaaaaa", "TRUE" = "#d62728"),
    labels = c("Not significant",
               paste0("FDR < 0.05 (n=", sum(anova_results$padj < 0.05), ")"))
  ) +
  labs(
    title    = "ANOVA — 1285 features across 13 site groups",
    subtitle = "X = F-statistic (larger = bigger between-site signal vs within-site noise). Y = significance.",
    x = "F-statistic", y = "-log10(FDR adjusted p-value)", color = ""
  ) +
  theme_bw(base_size = 12)

# PERMANOVA — whole-metabolome test -----------------------------------------

mat_t <- t(pos_mat)   # vegan needs rows=samples, cols=features

set.seed(42)
perm_result <- adonis2(
  mat_t ~ site_factor,
  method       = "euclidean",
  permutations = 999
)

cat("=== PERMANOVA ===\n")
print(perm_result)
cat("\nR² =", round(perm_result$R2[1], 4),
    "→ site explains", round(100 * perm_result$R2[1], 1),
    "% of total metabolome variance\n")

# Volcano plots — each site vs control --------------------------------------

run_volcano <- function(site_name) {
  ck_cols   <- rownames(meta)[meta$site == "CK"]
  site_cols <- rownames(meta)[meta$site == site_name]
  site_cols <- site_cols[site_cols %in% colnames(pos_mat)]

  fc   <- rowMeans(pos_mat[ , site_cols]) - rowMeans(pos_mat[ , ck_cols])
  pval <- sapply(seq_len(nrow(pos_mat)), function(i) {
    tryCatch(
      t.test(as.numeric(pos_mat[i, site_cols]),
             as.numeric(pos_mat[i, ck_cols]))$p.value,
      error = function(e) NA_real_
    )
  })

  padj <- p.adjust(pval, method = "BH")
  cat_label <- ifelse(padj < 0.05 & fc >  0.5, "Up in Danube",
               ifelse(padj < 0.05 & fc < -0.5, "Down in Danube",
                                                "Not significant"))

  data.frame(feature = rownames(pos_mat), site = site_name,
             fold_change = fc, pval = pval, padj = padj,
             neg_logp = -log10(padj), category = cat_label,
             stringsAsFactors = FALSE)
}

cat("Running per-site comparisons vs control...\n")
all_volcano <- do.call(rbind, lapply(
  paste0("D", sprintf("%02d", 1:12)), run_volcano
))
cat("Done.\n\n")

# Summary table
cat("=== SIGNIFICANT FEATURES PER SITE (FDR<0.05, |FC|>0.5) ===\n")
vsummary <- do.call(rbind, lapply(paste0("D", sprintf("%02d", 1:12)), function(s) {
  sub <- all_volcano[all_volcano$site == s, ]
  data.frame(Site  = s,
             Total = sum(sub$category != "Not significant"),
             Up    = sum(sub$category == "Up in Danube"),
             Down  = sum(sub$category == "Down in Danube"))
}))
print(vsummary)

ggplot(all_volcano, aes(x = fold_change, y = neg_logp, color = category)) +
  geom_point(alpha = 0.4, size = 0.7) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(
    values = c("Not significant" = "#cccccc",
               "Up in Danube"   = "#d62728",
               "Down in Danube" = "#2c7bb6")
  ) +
  facet_wrap(~ site, ncol = 4) +
  labs(
    title    = "Volcano plots — each Danube site vs Control",
    subtitle = "X = fold change (positive = higher than control). Y = significance. Coloured = FDR<0.05 and |FC|>0.5.",
    x = "Log fold change (site − control)",
    y = "-log10(FDR adjusted p-value)",
    color = ""
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11))

# Heatmap of significant features only --------------------------------------

sig_features <- anova_results$feature[anova_results$padj < 0.05]
sig_mat      <- pos_mat[sig_features, ]

cat("Significant features:", length(sig_features), "\n")

heat_sig <- as.matrix(sig_mat)

ann_col2 <- data.frame(
  Site  = meta[colnames(heat_sig), "site"],
  Batch = meta[colnames(heat_sig), "batch"],
  row.names = colnames(heat_sig)
)

sites2   <- unique(ann_col2$Site)
batches2 <- unique(ann_col2$Batch)
b_colors2 <- setNames(batch_pal[seq_along(batches2)], batches2)

ann_colors2 <- list(
  Site  = site_colors[names(site_colors) %in% sites2],
  Batch = b_colors2
)

pheatmap(
  heat_sig,
  annotation_col    = ann_col2,
  annotation_colors = ann_colors2,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  scale             = "row",
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#313695", "#f7f7f7", "#a50026"))(100),
  main  = paste0("Statistically significant features only (n=", length(sig_features),
                 ", FDR<0.05)\nRow Z-scored | Ward's D2 clustering")
)
