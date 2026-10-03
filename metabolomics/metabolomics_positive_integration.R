# ============================================================================
# POSITIVE MODE METABOLOMICS — PATHWAY ANNOTATION & TEAM HANDOFF
# Daphnia magna × Danube River multi-omics project
# ============================================================================
#
# Prerequisite: Run metabolomics_positive_EDA.R first (or load its workspace)
#
# Input files:
#   metabolomics_workspace.RData             — from EDA script
#   polar_pos_pkl_to_kegg_annotations.tsv    — KEGG compound annotations
#   polar_pos_pkl_to_hmdb_annotations.tsv    — HMDB compound annotations
#   anova_posthoc.csv                        — MetaboAnalyst ANOVA results
#
# Outputs:
#   sig_kegg_ids_site.txt                    — KEGG IDs from R ANOVA
#   sig_features_site_annotated.csv          — annotated significant features (R)
#   sig_features_pos_annotated.csv           — annotated significant features (MetaboAnalyst)
#   pos_sig_kegg_ids.txt                     — KEGG IDs for pathway analysis
#   pos_sig_hmdb_ids.txt                     — HMDB IDs for MSEA
#   pos_background_hmdb_ids.txt              — background HMDB IDs for MSEA
#   metabolomics_complete.RData              — full workspace for team
#   significant_features_annotated_pos.csv   — 100 sig features + annotations
#   site_significant_counts_pos.csv          — per-site summary
# ============================================================================

# Libraries -----------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

# Load workspace from EDA script -------------------------------------------

load("metabolomics_workspace.RData")

# NOTE: anova_results, sig_features, sig_mat, pca_df, pca_var, pca_res,
#       feat_var, feat_mean, all_volcano, vsummary are computed in the EDA
#       script. If running standalone, run the EDA script first or load
#       metabolomics_complete.RData (produced at the end of this script).

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

# ============================================================================
# ANNOTATION MATCHING — R ANOVA RESULTS
# Uses your own R ANOVA results, NOT MetaboAnalyst output
# MetaboAnalyst renames features internally so IDs don't match
# ============================================================================

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

# ============================================================================
# ANNOTATION MATCHING — MetaboAnalyst ANOVA results
# Matches on mz value (tolerance 0.001 Da) because MetaboAnalyst
# renames features internally so string IDs don't match
# ============================================================================

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

# ---- Match on mz value within tolerance 0.001 Da ----

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

# ---- KEGG IDs for pathway analysis ----

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

# ============================================================================
# HMDB IDs for MetaboAnalyst MSEA
# ============================================================================

sig_hmdb_expanded <- sig %>%
  filter(!is.na(hmdb_ids) & hmdb_ids != "") %>%
  separate_rows(hmdb_ids, sep = ";") %>%
  mutate(hmdb_ids = trimws(hmdb_ids)) %>%
  filter(nchar(hmdb_ids) > 0)

hmdb_list <- unique(sig_hmdb_expanded$hmdb_ids)
hmdb_list <- as.character(hmdb_list[nchar(hmdb_list) > 0])

cat("=== HMDB IDs for MetaboAnalyst MSEA ===\n")
cat("Total unique HMDB IDs:", length(hmdb_list), "\n\n")

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
print(summary_df)

# Session info --------------------------------------------------------------

sessionInfo()
