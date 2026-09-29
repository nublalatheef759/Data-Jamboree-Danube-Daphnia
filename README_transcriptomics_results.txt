SHARED TRANSCRIPTOMICS RESULTS
==============================

Summary plots and tables for everyone.

SITE RANKINGS:
- per_site_10x_summary.csv / per_site_1x_summary.csv
- degs_per_site.png : bar plot

KEY FINDINGS:
- D11 (1570 DEGs) and D12 (1605 DEGs) at 1x are the most affected sites
- D02, D03, D04, D09 dominate at 10x
- D08 at 10x: 0 DEGs (anomaly)

WGCNA MODULES:
- module_dendrogram.png : 17 modules from 5000 most variable genes
- module_chemical_heatmap.png : KEY FIGURE — module-chemical correlations
- module_trait_metadata.png : module vs site/REF
- meblue_eigengene_boxplot.png : MEblue (Diuron-linked) across sites

CHEMISTRY AT D11/D12:
- top_chems_per_site.png : wastewater chemicals concentrating downstream
- all_chemicals_d11_d12_summary.csv : full breakdown

DEG GENERALITY:
- deg_overlap_heatmap.png : Jaccard between site-dose conditions
- deg_frequency_histogram.png : universal vs specific DEGs

SUMMARY FOR SLIDES:
- 2 distinct pollution regimes along the Danube:
  (a) Agricultural/personal-care chemistry upstream-middle (Diuron, DEET)
  (b) Municipal wastewater downstream (Metformin, sweeteners, drugs)
- Wastewater signature at D11/D12 is associated with the largest biological response
- D12 also has a unique pesticide pulse (Chlorpropham, Carbaryl, Diazinon,
  Simazine, Chlorophene) — found ONLY at D12
- 29 universal pollution responders cluster into:
  * Lipid metabolism (NPC2, FABP12, LIPC)
  * Phase II detox (SULT1A5)
  * Digestive proteases (suppressed)
