MULTI-OMICS INTEGRATION INPUTS
==============================

For integrating transcriptomics with metabolomics (e.g. via MOFA, RGCCA, SGCCA).

PRIMARY INTEGRATION FEATURES:
- module_eigengenes.csv : 17 WGCNA module eigengenes per sample (148 samples)
  RECOMMENDED transcriptomic input for MOFA — reduces dimensionality
  from 5000 genes to ~17 modules, each with biological meaning.

GENE-LEVEL DATA (if needed):
- vst_top5000_for_iDEP.csv : top 5000 variable VST genes x 148 samples
  Use only if your method requires gene-level features.

CONTEXT FOR INTERPRETATION:
- module_gene_assignments.csv : which genes are in each module
- module_chemical_correlations.csv : how each module correlates with the 38
  chemicals — helps you link transcriptomic factors to chemistry
- sample_sheet.csv : metadata (148 rows: SampleID, REF, Site, Description)

SAMPLE EXCLUSIONS:
2 samples removed for low library size: D03A6, D03B3
Final n = 148 samples

KEY MODULES BY CHEMICAL ASSOCIATION:
- MEblue (741 genes)  : strongest correlation with Diuron (r=+0.585)
                        Also responds to wastewater chemicals (Acesulfame, Sucralose,
                        Benzotriazoles). Enriched for detox/lipid metabolism.
- MEpink, MEmagenta   : co-respond with MEblue, similar chemistry
- MEpurple            : negatively correlates with pollution (suppressed biology)

