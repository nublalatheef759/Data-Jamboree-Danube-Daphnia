Chemical Correlation — README

Overview
Spearman correlation analysis between hub gene expression (VST counts) and 90 water chemical concentrations across 12 Danube sites (D01–D12). This links the network analysis hub genes to specific pollutants driving their expression.

Inputs
- rna_vst_counts.csv — VST-normalised gene expression matrix (from Rayane)
- water_chemicals.tsv — 90 chemical concentrations measured at each site
- rna_dma_to_hsa_mappings.tsv — Daphnia magna to human ortholog mappings (to map hub gene symbols back to Daphnia gene IDs in the VST matrix)

Method
1. Mapped 9 hub genes (ITGB1, ANPEP, CD36, ACE, CD81, ABCA1, ABCG2, HSP90AA1, CTSD) back to Daphnia gene IDs using the ortholog mapping file
2. Extracted hub gene expression from VST matrix (averaged across Daphnia paralogs per hub gene)
3. Averaged expression per site (D01–D12), excluding control samples (no chemical data for controls)
4. Computed Spearman rank correlation between each hub gene and each chemical (9 genes x 90 chemicals = 810 tests)
5. Filtered to significant correlations (p < 0.05) and generated a heatmap

Tool: R (tidyverse, pheatmap, RColorBrewer). Script: hub_gene_chemical_correlation.R

Files in this folder
- hub_gene_chemical_correlation_heatmap.pdf — Correlation heatmap of hub genes vs significant chemicals. Blue = negative correlation, red = positive, * = p<0.05, ** = p<0.01
- hub_gene_chemical_correlations_full.csv — All 810 correlations with Spearman r and p-values

Key findings

72 significant correlations (p < 0.05) across 9 hub genes and 32 chemicals.

Strongest correlations:
- CD36 x Mecoprop: r = 0.830, p = 0.0008 (herbicide)
- ABCG2 x 2,4-Dinitrophenol: r = 0.832, p = 0.0008 (industrial pollutant)
- CD36 x Diuron: r = 0.776, p = 0.003 (herbicide — independently confirms Rayane's WGCNA finding of MEblue-Diuron r = 0.585)
- ANPEP x Carbamazepine: r = 0.755, p = 0.005 (anticonvulsant drug)
- ANPEP x Diuron: r = 0.738, p = 0.006

High-CV chemicals (identified from EDA as most site-variable):
- Terbuthylazine — significantly correlated with ABCG2 (r=0.682), ACE (r=0.654), ANPEP (r=0.650)
- Oxazepam — significantly correlated with CD36 (r=0.579)
- Sulfamethoxazole — no significant hub gene correlations
- Metoprolol — no significant hub gene correlations

Biological interpretation:
- CD36 (lipid receptor, #3 hub) is most strongly driven by herbicides (diuron, mecoprop, chlorotoluron, isoproturon)
- ABCG2 (xenobiotic efflux pump) responds to industrial chemicals (2,4-dinitrophenol) and the pesticide terbuthylazine
- ANPEP (#2 hub) responds broadly to pharmaceuticals (carbamazepine) and herbicides (diuron)
- Terbuthylazine drives 3 hub genes simultaneously, making it a candidate key stressor
- Diuron correlates with both CD36 and ANPEP, confirming its role as the strongest single pollution driver (consistent with Rayane's WGCNA)

