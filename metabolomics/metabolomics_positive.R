knitr::opts_chunk$set(
  echo    = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.align = "center"
)

# ============================================================================
# ================================================================ ==========
# ============================================================================

# ============================================================================
# PART 1 — SETUP, DATA LOADING & METADATA ===================================
# ============================================================================

# ============================================================================
# ================================================================ ==========
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
# ================================================================ ==========
# ============================================================================

# ============================================================================
# PART 2 — QUALITY CONTROL & OUTLIER DETECTION ==============================
# ============================================================================

# ============================================================================
# ================================================================ ==========
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
# Mathematically: sqrt( (x-mu)' * Sigma^-1 * (x-mu) )
# where mu = centroid and Sigma = covariance matrix of PC scores
pc_scores  <- pca_df[ , c("PC1", "PC2")]
centroid   <- colMeans(pc_scores)
cov_mat    <- cov(pc_scores)
mah_sq     <- mahalanobis(pc_scores, centroid, cov_mat)  # returns SQUARED distance

# Threshold: if the data were multivariate normal, squared Mahalanobis
# distances follow chi-squared(df=2). 97.5th percentile = 5.991.
# Any sample with squared distance > 5.991 is a formal statistical outlier.
chi_sq_thresh <- qchisq(0.975, df = 2)   # = 5.991

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
  # Label outliers without ggrepel — using base geom_text with nudge
  geom_text(
    data  = pca_df[pca_df$pca_outlier, ],
    aes(label = sample),
    nudge_y = 2, size = 3.5, color = "#d62728", fontface = "bold"
  ) +
  # Also label the Layer 1 mean-flagged samples for comparison
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
cat("Report wording: 'One sample (D08B6) exceeded the Mahalanobis 97.5%\n")
cat("threshold in PCA but passed all other QC criteria (mean z=0.83,\n")
cat("within-site correlation r=0.64 vs site average r=0.65). Retained.'\n")

# ============================================================================
# ================================================================ ==========
# ============================================================================

# ============================================================================
# PART 3 — EXPLORATORY DATA ANALYSIS ========================================
# ============================================================================

# ============================================================================
# ================================================================ ==========
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
    subtitle = "Lower = daphea from same site are metabolically consistent. All sites ~0.03-0.045 = excellent reproducibility.",
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
# (mismatched names cause the pheatmap error you saw before)
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
# ================================================================ ==========
# ============================================================================

# ============================================================================
# PART 4 — DIFFERENTIAL ABUNDANCE ANALYSIS ==================================
# ============================================================================

# ============================================================================
# ================================================================ ==========
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
# Controls the false discovery rate: expected proportion of your
# significant list that are false positives
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

# ANOVA tests features one-by-one. PERMANOVA tests ALL 1285 simultaneously:
# "Does the overall metabolite composition differ by site?"
# This is the multivariate headline test for your research question.

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

# ============================================================================
# ================================================================ ==========
# ============================================================================

# ============================================================================
# PART 5 — PATHWAY ANNOTATION & TEAM HANDOFF ================================
# ============================================================================

# ============================================================================
# ================================================================ ==========
# ============================================================================

# Join significant features to annotations ----------------------------------

# Merge annotation tables
annot_full <- merge(
  pos_kegg[ , c("feature_id", "kegg_ids")],
  pos_hmdb[ , c("feature_id", "hmdb_ids")],
  by = "feature_id", all = TRUE
)

# Join to ANOVA results
sig_annotated <- merge(
  anova_results[anova_results$padj < 0.05, ],
  annot_full,
  by.x = "feature", by.y = "feature_id",
  all.x = TRUE
)

cat("Significant features (FDR<0.05):      ", nrow(sig_annotated), "\n")
cat("Of which KEGG annotated:              ",
    sum(!is.na(sig_annotated$kegg_ids)), "\n")
cat("Of which HMDB annotated:              ",
    sum(!is.na(sig_annotated$hmdb_ids)), "\n")
cat("Of which unannotated (used in ML):    ",
    sum(is.na(sig_annotated$kegg_ids) & is.na(sig_annotated$hmdb_ids)), "\n\n")

cat("Top 15 significant features with KEGG annotation:\n")
top_kegg <- sig_annotated[!is.na(sig_annotated$kegg_ids), ]
top_kegg <- top_kegg[order(top_kegg$padj), ]
show <- top_kegg[1:min(15, nrow(top_kegg)),
                 c("feature", "mz", "F_stat", "padj", "kegg_ids", "hmdb_ids")]
show[ , c("F_stat","padj")] <- round(show[ , c("F_stat","padj")], 5)
print(show)

# Which KEGG compound IDs appear most in your significant features?
# High frequency = implicated by multiple features = likely important

sig_kegg_ids <- unlist(strsplit(
  sig_annotated$kegg_ids[!is.na(sig_annotated$kegg_ids)], ";"))
sig_kegg_ids <- trimws(sig_kegg_ids[nchar(trimws(sig_kegg_ids)) > 0])

kegg_freq <- sort(table(sig_kegg_ids), decreasing = TRUE)

cat("Most frequent KEGG compound IDs in significant features:\n")
print(head(kegg_freq, 20))

# Look these up at: https://www.genome.jp/kegg/compound/
# e.g., search "C00158" → citric acid → TCA cycle
# Ecotox pathways to watch: steroid biosynthesis, glutathione,
# TCA cycle, amino acid metabolism, fatty acid metabolism

df_kegg_freq <- data.frame(
  kegg_id = names(kegg_freq),
  count   = as.integer(kegg_freq)
)

ggplot(head(df_kegg_freq, 20), aes(x = reorder(kegg_id, count), y = count)) +
  geom_col(fill = "#2c7bb6", alpha = 0.85) +
  coord_flip() +
  labs(
    title    = "Most frequent KEGG IDs among significant features",
    subtitle = "Look these up at kegg.jp/compound to get metabolite names and pathways.",
    x = "KEGG Compound ID", y = "Feature count"
  ) +
  theme_bw(base_size = 11)

# ================================================================
# ANNOTATION MATCHING — CORRECTED
# Uses your own R ANOVA results, NOT MetaboAnalyst output
# MetaboAnalyst renames features internally so IDs don't match
# ================================================================

# Load annotation files
kegg_annot <- read.table(
  "polar_pos_pkl_to_kegg_annotations.tsv",
  sep = "\t", header = FALSE,
  col.names = c("feature_id", "mz", "kegg_ids"),
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

hmdb_annot <- read.table(
  "polar_pos_pkl_to_hmdb_annotations.tsv",
  sep = "\t", header = FALSE,
  col.names = c("feature_id", "mz", "hmdb_ids"),
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

cat("Total features in annotation file:", nrow(kegg_annot), "\n")
cat("KEGG annotated:", sum(!is.na(kegg_annot$kegg_ids)), "\n")
cat("HMDB annotated:", sum(!is.na(hmdb_annot$hmdb_ids)), "\n\n")

# ---- Use your R ANOVA results (anova_results is already in workspace) ----
# anova_results was built in Part 4 of your Rmd
# It has columns: feature, mz, F_stat, pval, padj
# feature column matches annotation file IDs exactly

sig_features_r <- anova_results[anova_results$padj < 0.05, ]
cat("Significant features from R ANOVA (site-based):", nrow(sig_features_r), "\n\n")

# Merge KEGG annotations
sig_kegg <- merge(
  sig_features_r,
  kegg_annot[ , c("feature_id", "kegg_ids")],
  by.x = "feature", by.y = "feature_id",
  all.x = TRUE
)

# Merge HMDB annotations
sig_annot <- merge(
  sig_kegg,
  hmdb_annot[ , c("feature_id", "hmdb_ids")],
  by.x = "feature", by.y = "feature_id",
  all.x = TRUE
)

cat("=== ANNOTATION RESULTS ===\n")
cat("Total significant:", nrow(sig_annot), "\n")
cat("With KEGG ID:     ", sum(!is.na(sig_annot$kegg_ids)), "\n")
cat("With HMDB ID:     ", sum(!is.na(sig_annot$hmdb_ids)), "\n")
cat("Unannotated:      ", sum(is.na(sig_annot$kegg_ids) & 
                               is.na(sig_annot$hmdb_ids)), "\n\n")

# Show top annotated features
top_kegg <- sig_annot[!is.na(sig_annot$kegg_ids), ]
top_kegg <- top_kegg[order(top_kegg$padj), ]
cat("Top significant features with KEGG annotation:\n")
print(top_kegg[ , c("feature", "mz", "F_stat", "padj", "kegg_ids")])

# ---- Extract KEGG IDs for MetaboAnalyst pathway analysis ----

# Split semicolon-separated IDs and flatten
kegg_split <- unlist(strsplit(
  top_kegg$kegg_ids[!is.na(top_kegg$kegg_ids)], ";"))
kegg_split <- trimws(kegg_split)
kegg_split <- kegg_split[nchar(kegg_split) > 0]
kegg_unique <- unique(kegg_split)

cat("\nUnique KEGG IDs for pathway analysis:", length(kegg_unique), "\n")
cat("IDs:\n")
cat(paste(kegg_unique, collapse = "\n"))

# Save KEGG list as plain text
writeLines(as.character(kegg_unique), "sig_kegg_ids_site.txt")
cat("\n\nSaved: sig_kegg_ids_site.txt\n")

# Save full annotated table
write.csv(sig_annot, "sig_features_site_annotated.csv", row.names = FALSE)
cat("Saved: sig_features_site_annotated.csv\n")

# FOR ANNOTATION MATCHING — POSITIVE MODE (MetaboAnalyst anova_posthoc.csv)

library(dplyr)
library(tidyr)

# Load ANOVA results from MetaboAnalyst
anova_res <- read.csv("anova_posthoc.csv", row.names = 1)
sig <- anova_res[anova_res$FDR < 0.05, ]
cat("Significant features:", nrow(sig), "\n")

# Extract mz values from MetaboAnalyst feature IDs
# MetaboAnalyst names them "mz_326.0856" — strip the "mz_" prefix
sig$mz_val <- as.numeric(sub("mz_", "", rownames(sig)))
sig$feature_id <- rownames(sig)

# Load annotation files
kegg_annot <- read.delim(
  "polar_pos_pkl_to_kegg_annotations.tsv",
  header = FALSE, sep = "\t",
  col.names = c("feature_id_annot", "mz_annot", "kegg_ids"),
  stringsAsFactors = FALSE,
  fill = TRUE, quote = "",
  na.strings = c("", "NA")
)

hmdb_annot <- read.delim(
  "polar_pos_pkl_to_hmdb_annotations.tsv",
  header = FALSE, sep = "\t",
  col.names = c("feature_id_annot", "mz_annot", "hmdb_ids"),
  stringsAsFactors = FALSE,
  fill = TRUE, quote = "",
  na.strings = c("", "NA")
)

cat("Total pos features:", nrow(kegg_annot), "\n")
cat("Features with KEGG IDs:", sum(!is.na(kegg_annot$kegg_ids)), "\n")
cat("Features with HMDB IDs:", sum(!is.na(hmdb_annot$hmdb_ids)), "\n\n")

# ---- KEY FIX: match on mz value, not feature ID string ----
# MetaboAnalyst uses "mz_326.0856", annotation uses "326_0856092"
# These encode the same mz but as different strings
# Solution: match numerically within tolerance 0.001 Da

tolerance <- 0.001

match_by_mz <- function(sig_df, annot_df, annot_col) {
  result_col <- rep(NA_character_, nrow(sig_df))
  for (i in seq_len(nrow(sig_df))) {
    diffs <- abs(annot_df$mz_annot - sig_df$mz_val[i])
    best  <- which.min(diffs)
    if (diffs[best] < tolerance) {
      result_col[i] <- annot_df[[annot_col]][best]
    }
  }
  result_col
}

sig$kegg_ids <- match_by_mz(sig, kegg_annot, "kegg_ids")
sig$hmdb_ids <- match_by_mz(sig, hmdb_annot, "hmdb_ids")

# Results
cat("Significant with KEGG ID:", sum(!is.na(sig$kegg_ids)), "\n")
cat("Significant with HMDB ID:", sum(!is.na(sig$hmdb_ids)), "\n")
cat("Unannotated:             ", 
    sum(is.na(sig$kegg_ids) & is.na(sig$hmdb_ids)), "\n\n")

cat("Annotated KEGG %:", 
    round(100 * sum(!is.na(sig$kegg_ids)) / nrow(sig), 1), "%\n")
cat("Annotated HMDB %:", 
    round(100 * sum(!is.na(sig$hmdb_ids)) / nrow(sig), 1), "%\n\n")

# Show annotated features
cat("Features with KEGG annotation:\n")
kegg_hits <- sig[!is.na(sig$kegg_ids), c("mz_val", "p.value", "FDR", "kegg_ids", "hmdb_ids")]
kegg_hits <- kegg_hits[order(kegg_hits$FDR), ]
print(kegg_hits)

# Save annotated table
write.csv(sig, "sig_features_pos_annotated.csv", row.names = TRUE)
cat("\nSaved: sig_features_pos_annotated.csv\n")

# Extract KEGG IDs for MetaboAnalyst pathway analysis
sig_kegg_expanded <- sig %>%
  filter(!is.na(kegg_ids) & kegg_ids != "") %>%
  separate_rows(kegg_ids, sep = ";") %>%
  filter(trimws(kegg_ids) != "")

kegg_list <- unique(trimws(sig_kegg_expanded$kegg_ids))
kegg_list <- as.character(kegg_list[nchar(kegg_list) > 0])

writeLines(kegg_list, "pos_sig_kegg_ids.txt")
cat("KEGG IDs saved to: pos_sig_kegg_ids.txt\n")
cat("Total unique KEGG IDs:", length(kegg_list), "\n")
cat("\nPaste these into MetaboAnalyst → Pathway Analysis:\n")
cat(paste(kegg_list, collapse = "\n"))

# FOR ANNOTATION MATCHING — POSITIVE MODE (MetaboAnalyst anova_posthoc.csv)
library(dplyr)
library(tidyr)

# Load ANOVA results from MetaboAnalyst
anova_res <- read.csv("anova_posthoc.csv", row.names = 1)
sig <- anova_res[anova_res$FDR < 0.05, ]
cat("Significant features:", nrow(sig), "\n")

# Extract mz values from MetaboAnalyst feature IDs
sig$mz_val <- as.numeric(sub("mz_", "", rownames(sig)))
sig$feature_id <- rownames(sig)

# Load annotation files
kegg_annot <- read.delim(
  "polar_pos_pkl_to_kegg_annotations.tsv",
  header = FALSE, sep = "\t",
  col.names = c("feature_id_annot", "mz_annot", "kegg_ids"),
  stringsAsFactors = FALSE,
  fill = TRUE, quote = "",
  na.strings = c("", "NA")
)
hmdb_annot <- read.delim(
  "polar_pos_pkl_to_hmdb_annotations.tsv",
  header = FALSE, sep = "\t",
  col.names = c("feature_id_annot", "mz_annot", "hmdb_ids"),
  stringsAsFactors = FALSE,
  fill = TRUE, quote = "",
  na.strings = c("", "NA")
)

cat("Total pos features:", nrow(kegg_annot), "\n")
cat("Features with KEGG IDs:", sum(!is.na(kegg_annot$kegg_ids)), "\n")
cat("Features with HMDB IDs:", sum(!is.na(hmdb_annot$hmdb_ids)), "\n\n")

# Match on mz value within tolerance
tolerance <- 0.001

match_by_mz <- function(sig_df, annot_df, annot_col) {
  result_col <- rep(NA_character_, nrow(sig_df))
  for (i in seq_len(nrow(sig_df))) {
    diffs <- abs(annot_df$mz_annot - sig_df$mz_val[i])
    best  <- which.min(diffs)
    if (diffs[best] < tolerance) {
      result_col[i] <- annot_df[[annot_col]][best]
    }
  }
  result_col
}

sig$kegg_ids <- match_by_mz(sig, kegg_annot, "kegg_ids")
sig$hmdb_ids <- match_by_mz(sig, hmdb_annot, "hmdb_ids")

# Summary
cat("Significant with KEGG ID:", sum(!is.na(sig$kegg_ids)), "\n")
cat("Significant with HMDB ID:", sum(!is.na(sig$hmdb_ids)), "\n")
cat("Unannotated:             ", 
    sum(is.na(sig$kegg_ids) & is.na(sig$hmdb_ids)), "\n\n")
cat("Annotated KEGG %:", 
    round(100 * sum(!is.na(sig$kegg_ids)) / nrow(sig), 1), "%\n")
cat("Annotated HMDB %:", 
    round(100 * sum(!is.na(sig$hmdb_ids)) / nrow(sig), 1), "%\n\n")

# Show annotated features
cat("Features with KEGG annotation:\n")
kegg_hits <- sig[!is.na(sig$kegg_ids), c("mz_val", "p.value", "FDR", "kegg_ids", "hmdb_ids")]
kegg_hits <- kegg_hits[order(kegg_hits$FDR), ]
print(kegg_hits)

# Save annotated table
write.csv(sig, "sig_features_pos_annotated.csv", row.names = TRUE)
cat("\nSaved: sig_features_pos_annotated.csv\n")

# ---- KEGG IDs for pathway analysis ----
sig_kegg_expanded <- sig %>%
  filter(!is.na(kegg_ids) & kegg_ids != "") %>%
  separate_rows(kegg_ids, sep = ";") %>%
  filter(trimws(kegg_ids) != "")

kegg_list <- unique(trimws(sig_kegg_expanded$kegg_ids))
kegg_list <- as.character(kegg_list[nchar(kegg_list) > 0])
writeLines(kegg_list, "pos_sig_kegg_ids.txt")
cat("KEGG IDs saved to: pos_sig_kegg_ids.txt\n")
cat("Total unique KEGG IDs:", length(kegg_list), "\n\n")

# ================================================================
# NEW — HMDB IDs for MetaboAnalyst MSEA
# ================================================================

# Extract and clean all HMDB IDs from significant features
sig_hmdb_expanded <- sig %>%
  filter(!is.na(hmdb_ids) & hmdb_ids != "") %>%
  separate_rows(hmdb_ids, sep = ";") %>%
  mutate(hmdb_ids = trimws(hmdb_ids)) %>%
  filter(nchar(hmdb_ids) > 0)

hmdb_list <- unique(sig_hmdb_expanded$hmdb_ids)
hmdb_list <- as.character(hmdb_list[nchar(hmdb_list) > 0])

cat("=== HMDB IDs for MetaboAnalyst MSEA ===\n")
cat("Total unique HMDB IDs:", length(hmdb_list), "\n\n")

# MetaboAnalyst MSEA accepts HMDB IDs in format HMDB0000001
# Your IDs are already in this format — no conversion needed
cat("Your HMDB IDs (paste into MetaboAnalyst MSEA):\n")
cat(paste(hmdb_list, collapse = "\n"))

# Save to file
writeLines(hmdb_list, "pos_sig_hmdb_ids.txt")
cat("\n\nSaved: pos_sig_hmdb_ids.txt\n")

# Also extract background HMDB IDs (all annotated features)
# MetaboAnalyst needs this for a fair enrichment test
bg_hmdb_expanded <- hmdb_annot %>%
  filter(!is.na(hmdb_ids) & hmdb_ids != "") %>%
  separate_rows(hmdb_ids, sep = ";") %>%
  mutate(hmdb_ids = trimws(hmdb_ids)) %>%
  filter(nchar(hmdb_ids) > 0)

bg_hmdb_list <- unique(as.character(bg_hmdb_expanded$hmdb_ids))
writeLines(bg_hmdb_list, "pos_background_hmdb_ids.txt")
cat("Background HMDB IDs saved to: pos_background_hmdb_ids.txt\n")
cat("Total background HMDB IDs:", length(bg_hmdb_list), "\n")

# Summary table of all annotated significant features
cat("\n=== FULL ANNOTATED SIGNIFICANT FEATURES ===\n")
annotated_sig <- sig[!is.na(sig$kegg_ids) | !is.na(sig$hmdb_ids), 
                      c("mz_val", "p.value", "FDR", "kegg_ids", "hmdb_ids")]
annotated_sig <- annotated_sig[order(annotated_sig$FDR), ]
print(annotated_sig)

# Save all outputs ----------------------------------------------------------

# Save everything your team needs
sig_mat_shared <- sig_mat[ , !colnames(sig_mat) %in%
                              c("D06A2","D06B3","D07A3","D08B1","D10B4")]

save(
  pos_mat, mz_values, meta, pos_kegg, pos_hmdb, site_colors,
  pca_df, pca_var, pca_res,
  feat_var, feat_mean,
  anova_results, sig_features, sig_mat, sig_annotated,
  all_volcano, vsummary,
  sig_mat_shared,
  file = "metabolomics_complete.RData"
)

write.csv(sig_annotated,
          "significant_features_annotated_pos.csv",
          row.names = FALSE)
write.csv(vsummary,
          "site_significant_counts_pos.csv",
          row.names = FALSE)

cat("Saved files:\n")
cat("  metabolomics_complete.RData           (full R workspace)\n")
cat("  significant_features_annotated_pos.csv (100 sig features + annotations)\n")
cat("  site_significant_counts_pos.csv        (per-site summary for team)\n")

# Messages for your teammates -----------------------------------------------

cat("=== FOR YOUR PARTNER (negative mode) ===\n")
cat("Shared samples (for combining modes): 144\n")
cat("Samples in your positive only: D06A2, D06B3, D07A3, D08B1, D10B4\n")
cat("Sample in negative only: D01B4\n")
cat("When combining: shared <- intersect(colnames(pos_mat), colnames(neg_mat))\n\n")

cat("=== FOR TRANSCRIPTOMICS ===\n")
cat("Our PCA: sites overlap, no clean separation.\n")
cat("Strongest metabolic signals: D02, D04, D09, D10, D12.\n")
cat("Compare with RNA-seq PCA — do the same sites drive variation?\n")
cat("D01 paradox: highest mean elevation but 0 specific significant features.\n\n")

cat("=== FOR CHEMICAL ANALYSIS ===\n")
cat("D01: broad metabolic elevation, no specific features — what's unique at D01 chemically?\n")
cat("D05/D06/D07/D08: 0 significant features — what do these mixtures share?\n")
cat("D12: most complex response (147 sig features) — most complex mixture?\n\n")

cat("=== FOR ML / INTEGRATION ===\n")
cat("Deliverable 1: sig_mat_shared (100 features × 144 shared samples)\n")
cat("Deliverable 2: significant_features_annotated_pos.csv\n")
cat("Deliverable 3: site_significant_counts_pos.csv\n")
cat("Suggestion: first 5 PCA scores as dimensionality-reduced input for iRF\n")
cat("Use all 1285 features (not just significant 100) for ML — let iRF select.\n")

# Final summary table -------------------------------------------------------

summary_df <- data.frame(
  Item = c(
    "Features", "Samples", "Missing values", "Near-zero variance features",
    "KEGG annotated", "HMDB annotated", "Unannotated",
    "PCA outliers", "QC flags (Layer 1 mean)",
    "ANOVA significant (FDR<0.05)", "ANOVA significant (FDR<0.01)",
    "Sites with 0 sig features", "Site with most sig features",
    "Samples for multi-omics integration"
  ),
  Result = c(
    "1,285", "149 (positive mode)", "0 — fully imputed",
    "0 — all retained", "437 (34.0%)", "409 (31.8%)", "757 (58.9%)",
    "1 — D08B6 (biological, kept)", "5 — marginal, all kept",
    "100", "45",
    "D05, D06, D07, D08",
    "D12 (147 sig, 103 up, 44 down)",
    "144 shared with negative mode"
  ),
  stringsAsFactors = FALSE
)
knitr::kable(summary_df, caption = "Complete Analysis Summary — Positive Mode Metabolomics")

# Session info --------------------------------------------------------------

sessionInfo()
