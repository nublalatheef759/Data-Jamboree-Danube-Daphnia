knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

# ============================================================================
# 1. Differential expression analysis (DESeq2) ==============================
# ============================================================================

# 1.1 Setup and data loading ------------------------------------------------

library(DESeq2)

counts <- read.csv("rna_raw_counts.csv", row.names = 1, check.names = FALSE)
counts <- round(counts)
ss <- read.csv("sample_sheet.csv", stringsAsFactors = FALSE)

# Remove failed samples (library sizes <1M reads)
ss <- ss[!ss$SampleID %in% c("D03A6", "D03B3"), ]
counts <- counts[, colnames(counts) %in% ss$SampleID]
ss <- ss[match(colnames(counts), ss$SampleID), ]

# Create group factor (site × concentration)
ss$group <- factor(paste0(ss$Site, "_", ss$REF))

# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = ss, design = ~group)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

print(levels(dds$group))

# 1.2 Extract 10x vs Control for all sites ----------------------------------

sites <- paste0("D", sprintf("%02d", 1:12))

for (site in sites) {
  group_name <- paste0(site, "_10x")
  if (group_name %in% levels(dds$group)) {
    res <- results(dds, contrast = c("group", group_name, "Control_Control"))
    write.csv(as.data.frame(res), paste0(site, "_10x_vs_Control.csv"))
    sig_count <- sum(!is.na(res$padj) & res$padj < 0.05, na.rm = TRUE)
    cat(site, "10x:", sig_count, "DEGs (padj < 0.05)\n")
  }
}

# 1.3 Extract 1x vs Control for all sites -----------------------------------

for (site in sites) {
  group_name <- paste0(site, "_1x")
  if (group_name %in% levels(dds$group)) {
    res <- results(dds, contrast = c("group", group_name, "Control_Control"))
    write.csv(as.data.frame(res), paste0(site, "_1x_vs_Control.csv"))
    sig_count <- sum(!is.na(res$padj) & res$padj < 0.05, na.rm = TRUE)
    cat(site, "1x:", sig_count, "DEGs (padj < 0.05)\n")
  }
}

# ============================================================================
# 2. Ortholog mapping =======================================================
# ============================================================================

# 2.1 Map 10x DEGs to human, Drosophila, and D. pulex orthologs -------------

ortho <- read.csv("ortholog_mapping_anchored.csv")

for (site in sites) {
  file <- paste0(site, "_10x_vs_Control.csv")
  if (file.exists(file)) {
    degs <- read.csv(file, row.names = 1)
    sig <- degs[!is.na(degs$padj) & degs$padj < 0.05, ]
    sig$gene_id <- rownames(sig)
    merged <- merge(sig, ortho, by.x = "gene_id", by.y = "dma", all.x = TRUE)
    write.csv(merged, paste0(site, "_10x_orthologs.csv"), row.names = FALSE)

    human_ids <- unique(merged$hsa_ensembl[!is.na(merged$hsa_ensembl)])
    writeLines(human_ids, paste0(site, "_10x_human.txt"))

    fly_ids <- unique(merged$dme_flybase[!is.na(merged$dme_flybase)])
    writeLines(fly_ids, paste0(site, "_10x_fly.txt"))

    dpx_ids <- unique(merged$dpx_id[!is.na(merged$dpx_id)])
    writeLines(dpx_ids, paste0(site, "_10x_dpx.txt"))

    cat(site, ":", nrow(sig), "DEGs ->", length(human_ids), "human,",
        length(fly_ids), "fly,", length(dpx_ids), "dpx\n")
  }
}

# Combined file
all_orthologs <- data.frame()
for (site in sites) {
  file <- paste0(site, "_10x_orthologs.csv")
  if (file.exists(file)) {
    df <- read.csv(file)
    if (nrow(df) > 0) {
      df$site <- site
      all_orthologs <- rbind(all_orthologs, df)
    }
  }
}
write.csv(all_orthologs, "all_sites_10x_orthologs_combined.csv", row.names = FALSE)
cat("Combined file:", nrow(all_orthologs), "rows across", length(unique(all_orthologs$site)), "sites\n")

# 2.2 Map 1x DEGs to orthologs ----------------------------------------------

for (site in sites) {
  file <- paste0(site, "_1x_vs_Control.csv")
  if (file.exists(file)) {
    degs <- read.csv(file, row.names = 1)
    sig <- degs[!is.na(degs$padj) & degs$padj < 0.05, ]
    sig$gene_id <- rownames(sig)
    merged <- merge(sig, ortho, by.x = "gene_id", by.y = "dma", all.x = TRUE)
    write.csv(merged, paste0(site, "_1x_orthologs.csv"), row.names = FALSE)

    human_ids <- unique(merged$hsa_ensembl[!is.na(merged$hsa_ensembl)])
    writeLines(human_ids, paste0(site, "_1x_human.txt"))

    fly_ids <- unique(merged$dme_flybase[!is.na(merged$dme_flybase)])
    writeLines(fly_ids, paste0(site, "_1x_fly.txt"))

    dpx_ids <- unique(merged$dpx_id[!is.na(merged$dpx_id)])
    writeLines(dpx_ids, paste0(site, "_1x_dpx.txt"))

    cat(site, ":", nrow(sig), "DEGs ->", length(human_ids), "human,",
        length(fly_ids), "fly,", length(dpx_ids), "dpx\n")
  }
}

all_orthologs_1x <- data.frame()
for (site in sites) {
  file <- paste0(site, "_1x_orthologs.csv")
  if (file.exists(file)) {
    df <- read.csv(file)
    if (nrow(df) > 0) {
      df$site <- site
      all_orthologs_1x <- rbind(all_orthologs_1x, df)
    }
  }
}
write.csv(all_orthologs_1x, "all_sites_1x_orthologs_combined.csv", row.names = FALSE)
cat("Combined file:", nrow(all_orthologs_1x), "rows across", length(unique(all_orthologs_1x$site)), "sites\n")

# ============================================================================
# 3. Prepare data for iDEP and pathway enrichment ===========================
# ============================================================================

# 3.1 Convert VST matrix to human Ensembl IDs for iDEP ----------------------

library(matrixStats)

vst <- read.csv("rna_vst_counts.csv", row.names = 1, check.names = FALSE)
vst <- vst[, !colnames(vst) %in% c("D03A6", "D03B3")]

vst$dma <- rownames(vst)
merged <- merge(vst, ortho[, c("dma", "hsa_ensembl")], by = "dma")
merged <- merged[!is.na(merged$hsa_ensembl), ]

# If multiple Daphnia genes map to same human gene, keep highest variance
sample_cols <- colnames(merged)[!colnames(merged) %in% c("dma", "hsa_ensembl")]
merged$var <- rowVars(as.matrix(merged[, sample_cols]))
merged <- merged[order(-merged$var), ]
merged <- merged[!duplicated(merged$hsa_ensembl), ]

rownames(merged) <- merged$hsa_ensembl
vst_human <- merged[, sample_cols]

write.csv(vst_human, "vst_human_for_idep.csv")
cat("Matrix:", nrow(vst_human), "genes x", ncol(vst_human), "samples\n")

# 3.2 Create iDEP design file -----------------------------------------------

ss <- read.csv("sample_sheet.csv", stringsAsFactors = FALSE)
ss <- ss[!ss$SampleID %in% c("D03A6", "D03B3"), ]

vst <- read.csv("vst_human_for_idep.csv", row.names = 1, check.names = FALSE)
ss <- ss[ss$SampleID %in% colnames(vst), ]

idep_design <- data.frame(matrix(ncol = nrow(ss), nrow = 1))
colnames(idep_design) <- ss$SampleID
idep_design[1, ] <- ss$REF
rownames(idep_design) <- "Concentration"

write.csv(idep_design, "idep_design.csv")
cat("Design file:", ncol(idep_design), "samples\n")

# 3.3 Drosophila FlyBase to Entrez ID conversion (for WebGestalt) -----------

library(org.Dm.eg.db)

fly <- readLines("D11_1x_fly.txt")
mapping <- mapIds(org.Dm.eg.db, keys = fly, keytype = "FLYBASE", column = "ENTREZID")
entrez_ids <- unique(na.omit(mapping))
writeLines(entrez_ids, "D11_1x_fly_entrez.txt")
cat("Mapped:", length(entrez_ids), "out of", length(fly), "\n")

# 3.4 Prepare pooled DEG lists for iDEP -------------------------------------

# Pooled 10x
degs <- read.csv("A_pooled_10x_vs_Control.csv")
merged <- merge(degs, ortho, by.x = "gene_id", by.y = "dma")
merged <- merged[!is.na(merged$hsa_ensembl), ]
merged <- merged[!duplicated(merged$hsa_ensembl), ]

idep_input <- data.frame(
  gene = merged$hsa_ensembl,
  log2FoldChange = merged$log2FoldChange,
  padj = merged$padj
)
write.csv(idep_input, "pooled_10x_for_idep.csv", row.names = FALSE)
cat("Pooled 10x genes:", nrow(idep_input), "\n")

# Pooled 1x
degs_1x <- read.csv("B_pooled_1x_vs_Control.csv")
merged_1x <- merge(degs_1x, ortho, by.x = "gene_id", by.y = "dma")
merged_1x <- merged_1x[!is.na(merged_1x$hsa_ensembl), ]
merged_1x <- merged_1x[!duplicated(merged_1x$hsa_ensembl), ]

idep_input_1x <- data.frame(
  gene = merged_1x$hsa_ensembl,
  log2FoldChange = merged_1x$log2FoldChange,
  padj = merged_1x$padj
)
write.csv(idep_input_1x, "pooled_1x_for_idep.csv", row.names = FALSE)
cat("Pooled 1x genes:", nrow(idep_input_1x), "\n")

# ============================================================================
# 4. STRING network input — MEblue module gene list =========================
# ============================================================================

meblue <- read.csv("meblue_genes_annotated.csv")
symbols <- meblue$hsa_symbol[!is.na(meblue$hsa_symbol) & meblue$hsa_symbol != "NA"]
writeLines(symbols, "meblue_human_symbols.txt")
cat("Genes for STRING:", length(symbols), "\n")

# ============================================================================
# 5. AOP diagram construction ===============================================
# ============================================================================

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

graph <- grViz("
digraph AOP {
  rankdir=TB
  graph [bgcolor='white', pad=0.5]
  node [shape=box, style='filled,rounded', fontname='Arial',
        fontsize=11, margin=0.3, penwidth=0.5]
  edge [color='#555555', penwidth=1.2, arrowsize=0.8]

  S   [label='Stressor\\nDanube pollution mixture\\nDiuron, metformin, sweeteners, pesticides', fillcolor='#F09595']
  MIE [label='MIE: Membrane receptor disruption\\nITGB1 (#1 hub), CD36 (#3 hub), CD81', fillcolor='#F0997B']
  KE1 [label='KE1: Stress response & xenobiotic efflux\\nHSP90AA1, ACE, SOD3, CAT\\nABCG2, ABCB1, ABCA1', fillcolor='#FAC775']
  KE2 [label='KE2: Lipid & lysosomal disruption\\nNPC2 (universal), ANPEP (#2 hub)\\nCTSD, SMPD1, GLB1', fillcolor='#9FE1CB']
  KE3 [label='KE3: Suppressed growth & reproduction\\n~10:1 down:up gene regulation\\nD11: 1570 DEGs | D12: 1605 DEGs', fillcolor='#CECBF6']
  AO  [label='Adverse Outcome\\nReduced fitness, fecundity & survival\\nDaphnia magna population decline', fillcolor='#D3D1C7']

  S -> MIE -> KE1 -> KE2 -> KE3 -> AO
}
")

svg_text <- export_svg(graph)
writeLines(svg_text, "AOP_diagram.svg")

# ============================================================================
# 6. Hub gene × chemical correlation analysis ===============================
# ============================================================================

# 6.1 Load data and define hub genes ----------------------------------------

library(tidyverse)
library(pheatmap)
library(RColorBrewer)

vst <- read.csv("rna_vst_counts.csv", row.names = 1, check.names = FALSE)
chem <- read.delim("water_chemicals.tsv", check.names = FALSE)
mapping <- read.delim("rna_dma_to_hsa_mappings.tsv", header = FALSE,
                       col.names = c("dma_id", "ensembl_ids"))

hub_genes <- data.frame(
  symbol = c("ITGB1", "ANPEP", "CD36", "ACE", "CD81",
             "ABCA1", "ABCG2", "HSP90AA1", "CTSD"),
  ensembl = c("ENSG00000150093", "ENSG00000166825", "ENSG00000135218",
              "ENSG00000159640", "ENSG00000110651", "ENSG00000165029",
              "ENSG00000118777", "ENSG00000080824", "ENSG00000117984"),
  methods = c("3/3", "3/3", "3/3", "3/3", "2/3",
              "2/3", "2/3", "2/3", "2/3")
)

# 6.2 Map hub genes to Daphnia IDs ------------------------------------------

get_dma_ids <- function(ensembl_id, mapping_df, vst_genes) {
  hits <- mapping_df %>%
    filter(str_detect(ensembl_ids, fixed(ensembl_id))) %>%
    pull(dma_id)
  hits[hits %in% vst_genes]
}

vst_genes <- rownames(vst)

hub_dma_map <- list()
for (i in 1:nrow(hub_genes)) {
  ids <- get_dma_ids(hub_genes$ensembl[i], mapping, vst_genes)
  hub_dma_map[[hub_genes$symbol[i]]] <- ids
  cat(hub_genes$symbol[i], ":", length(ids), "Daphnia orthologs\n")
}

# 6.3 Average expression per site -------------------------------------------

sample_names <- colnames(vst)
sample_sites <- data.frame(
  sample = sample_names,
  site = sub("([A-Z]\\d{2}).*", "\\1", sample_names)
)
sample_sites <- sample_sites[grepl("^D\\d{2}", sample_sites$site), ]

hub_site_expr <- matrix(NA, nrow = length(hub_dma_map), ncol = 12,
                         dimnames = list(names(hub_dma_map),
                                         paste0("D", sprintf("%02d", 1:12))))

for (gene_name in names(hub_dma_map)) {
  dma_ids <- hub_dma_map[[gene_name]]
  if (length(dma_ids) == 0) next
  gene_expr <- colMeans(vst[dma_ids, , drop = FALSE])
  for (site in colnames(hub_site_expr)) {
    site_samples <- sample_sites$sample[sample_sites$site == site]
    if (length(site_samples) > 0) {
      hub_site_expr[gene_name, site] <- mean(gene_expr[site_samples], na.rm = TRUE)
    }
  }
}

cat("\nHub gene expression matrix (sites):\n")
print(round(hub_site_expr, 2))

# 6.4 Compute Spearman correlations -----------------------------------------

chem_mat <- as.matrix(chem[, paste0("D", sprintf("%02d", 1:12))])
rownames(chem_mat) <- chem$ChemName
chem_mat[is.na(chem_mat)] <- 0

cor_mat <- matrix(NA, nrow = nrow(hub_site_expr), ncol = nrow(chem_mat),
                   dimnames = list(rownames(hub_site_expr), rownames(chem_mat)))
pval_mat <- cor_mat

for (i in 1:nrow(hub_site_expr)) {
  for (j in 1:nrow(chem_mat)) {
    test <- cor.test(hub_site_expr[i, ], chem_mat[j, ],
                      method = "spearman", exact = FALSE)
    cor_mat[i, j] <- test$estimate
    pval_mat[i, j] <- test$p.value
  }
}

# 6.5 Filter and generate heatmap -------------------------------------------

sig_chems <- apply(pval_mat, 2, function(x) any(x < 0.05, na.rm = TRUE))
cat("Chemicals with ≥1 significant hub gene correlation:", sum(sig_chems), "\n")

high_cv_chems <- c("Terbuthylazine", "Sulfamethoxazole", "Metoprolol", "Oxazepam")
high_cv_present <- high_cv_chems[high_cv_chems %in% rownames(chem_mat)]

selected_chems <- unique(c(rownames(chem_mat)[sig_chems], high_cv_present))

if (length(selected_chems) > 30) {
  min_pvals <- apply(pval_mat[, selected_chems], 2, min, na.rm = TRUE)
  selected_chems <- names(sort(min_pvals))[1:30]
  selected_chems <- unique(c(selected_chems, high_cv_present))
}

cor_subset <- cor_mat[, selected_chems, drop = FALSE]
pval_subset <- pval_mat[, selected_chems, drop = FALSE]

sig_stars <- matrix("", nrow = nrow(pval_subset), ncol = ncol(pval_subset))
sig_stars[pval_subset < 0.05] <- "*"
sig_stars[pval_subset < 0.01] <- "**"

pheatmap(t(cor_subset),
         display_numbers = t(sig_stars),
         color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         fontsize = 9,
         fontsize_number = 12,
         angle_col = 0,
         main = "Hub gene × chemical correlation (Spearman)\n* p<0.05, ** p<0.01")

# 6.6 High-CV chemical correlations -----------------------------------------

if (length(high_cv_present) > 0) {
  for (ch in high_cv_present) {
    cat("\n", ch, ":\n")
    for (g in rownames(cor_mat)) {
      cat(sprintf("  %s: r=%.3f, p=%.4f %s\n",
                  g, cor_mat[g, ch], pval_mat[g, ch],
                  ifelse(pval_mat[g, ch] < 0.05, "*", "")))
    }
  }
}

# 6.7 Export full results ---------------------------------------------------

results <- data.frame()
for (i in 1:nrow(cor_mat)) {
  for (j in 1:ncol(cor_mat)) {
    results <- rbind(results, data.frame(
      hub_gene = rownames(cor_mat)[i],
      chemical = colnames(cor_mat)[j],
      spearman_r = cor_mat[i, j],
      p_value = pval_mat[i, j]
    ))
  }
}

results <- results %>%
  mutate(significant = p_value < 0.05) %>%
  arrange(p_value)

write.csv(results, "hub_gene_chemical_correlations_full.csv", row.names = FALSE)

cat("\nTop 20 strongest hub gene-chemical correlations:\n")
print(head(results, 20))
