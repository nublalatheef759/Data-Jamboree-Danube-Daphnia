# Data Jamboree — Multi-Omics Integration (Corrected)
## Role: Multi-Omics Integration (MOFA2)
## Reanalysis by: [Your Name]

---

## Summary

MOFA2 was rerun after identifying that the original analysis used the wrong positive-mode metabolomics file (`polar_pos_pqn.csv` — PQN only, no glog transformation). The corrected analysis uses `polar_pos_pqn_imputed_glog.csv`, matching the preprocessing applied to the negative-mode data. This correction revealed shared variance across all three omics layers that was absent in the original model.

---

## Files in this folder

| File | Description |
|------|-------------|
| vst_top5000_for_iDEP.csv | Top 5000 most variable VST-normalised genes × 142 shared samples (from Rayane) |
| polar_neg_pqn_imputed_glog.csv | Negative mode metabolomics — PQN, KNN imputed, glog-transformed, 2331 features |
| polar_pos_pqn_imputed_glog.csv | Positive mode metabolomics — PQN, KNN imputed, glog-transformed, 1285 features |
| module_eigengenes.csv | 17 WGCNA module eigengenes × 148 samples (from Rayane) |
| ortholog_mapping_anchored.csv | 14,014 Dma gene ortholog lookup table (Dma → Dpx → Dme → Hsa) |
| sample_sheet.csv | Sample metadata — SampleID, REF, Site, Description |
| MOFA2_integration.R | Full R script used to run the analysis |
| mofa_model_corrected.hdf5 | Trained MOFA2 model (10 factors, 142 samples, 3 views) |
| mofa_variance_explained.svg | Heatmap of R² variance explained per factor per view |
| mofa_variance_total.svg | Total variance explained per view (bar chart) |
| mofa_scatter_concentration.svg | Factor 1 vs Factor 2 scatter coloured by concentration (Control/1x/10x) |
| mofa_scatter_site.svg | Factor 1 vs Factor 2 scatter coloured by site |
| mofa_f1_vs_f5.svg | Factor 1 vs Factor 5 scatter — shared multi-omics factor |
| mofa_factor5_by_site.svg | Boxplot of Factor 5 scores across Danube sites D01–D12 |
| mofa_wgcna_integration.svg | Pearson correlation heatmap — MOFA2 factors vs WGCNA module eigengenes |
| mofa_wgcna_integration.png | Same heatmap as PNG |
| mofa_top_rna_F1.svg | Top 15 RNA gene loadings on Factor 1 |
| mofa_top_rna_F5.svg | Top 15 RNA gene loadings on Factor 5 |
| mofa_top_neg_F1.svg | Top 15 negative metabolomics loadings on Factor 1 |
| mofa_top_pos_F1.svg | Top 15 positive metabolomics loadings on Factor 1 |
| top50_rna_factor1_annotated.csv | Top 50 Factor 1 RNA genes with human ortholog symbols |
| top50_rna_factor5_annotated.csv | Top 50 Factor 5 RNA genes with human ortholog symbols |
| MOFA_Slides_v3.pptx | Presentation slides (2 slides) |

---

## Methods

### 1. Sample Matching
- RNA (VST): 148 samples
- Negative metabolomics: 146 samples
- Positive metabolomics: 150 samples
- Intersection across all three layers = **142 shared samples**

### 2. MOFA2
- Platform: BlueBear HPC, University of Birmingham
- Software: MOFA2 v1.10.0 (R 4.3.1), mofapy2 backend via basilisk
- Views:
  - RNA: 5,000 most variable VST genes (not eigengenes)
  - Metabo_Neg: 2,331 glog-transformed features
  - Metabo_Pos: 1,285 glog-transformed features
- All views scaled to unit variance (`scale_views = TRUE`)
- Factors: 10, seed = 42, convergence mode = medium

### 3. Key Difference from Original Analysis
- Original used `polar_pos_pqn.csv` (PQN only, no imputation, no glog)
- Corrected uses `polar_pos_pqn_imputed_glog.csv` (PQN + KNN imputed + glog)
- Result: positive-mode metabolomics now contributes meaningful variance (up to 21.8% on Factor 3, vs max 3.33% in original)

---

## Results — Variance Explained per Factor

| Factor | RNA (%) | Metabo_Neg (%) | Metabo_Pos (%) | Interpretation |
|--------|---------|----------------|----------------|----------------|
| Factor 1 | 0.01 | 21.95 | 13.35 | Shared metabolomics — chemical uptake |
| Factor 2 | 26.15 | 0.16 | 0.12 | Transcriptomics — biological reprogramming |
| Factor 3 | 0.01 | 2.17 | 21.81 | Positive-mode metabolomics dominant |
| Factor 4 | 18.78 | 0.53 | 0.63 | RNA-specific |
| **Factor 5** | **7.53** | **4.59** | **5.08** | **Shared across all 3 layers** |
| Factor 6 | 0.01 | 14.25 | 0.08 | Negative-mode specific |
| Factor 7 | 0.03 | 9.47 | 1.07 | Negative-mode dominant |
| Factor 8 | 0.36 | 5.93 | 2.11 | Mixed metabolomics |
| Factor 9 | 6.12 | 0.01 | 0.02 | RNA-specific |
| Factor 10 | 3.32 | 0.05 | 0.31 | RNA-specific |

---

## Results — WGCNA Cross-Validation

| MOFA Factor | WGCNA Module | Pearson r | Interpretation |
|-------------|-------------|-----------|----------------|
| Factor 2 | MEturquoise | 0.99 | Largest transcriptomic module |
| **Factor 5** | **MEblue** | **0.98** | **Pollution-responsive module (Diuron-correlated)** |
| Factor 4 | MEblack | -0.88 | Inverse correlation |
| Factor 4 | MEbrown | -0.95 | Inverse correlation |

---

## Results — Top Factor 5 Genes (Human Orthologs)

| Daphnia Gene | Human Symbol | Function |
|-------------|-------------|----------|
| Dapma7bEVm015017 | CHI3L1 | Chitinase — chitin/exoskeleton metabolism |
| Dapma7bEVm010425 | AZU1 | Azurocidin — antimicrobial immune defence |
| Dapma7bEVm001244 | CTRB1 | Chymotrypsinogen — digestive enzyme |
| Dapma7bEVm000253 | CPA1 | Carboxypeptidase — digestive enzyme |
| Dapma7bEVm015167 | GM2A | Ganglioside metabolism |
| Dapma7bEVm028529 | NPC2 | Cholesterol transport |
| Dapma7bEVm028424 | PRSS45 | Serine protease |

---

## Key Findings

1. **RNA and metabolomics capture complementary biology** — Factor 1 (metabolomics) and Factor 2 (transcriptomics) are largely orthogonal
2. **Factor 5 is the shared multi-omics factor** — 7.5% RNA + 4.6% neg + 5.1% pos variance, the only factor with meaningful contribution from all three layers
3. **Factor 5 cross-validates with WGCNA MEblue** (r = 0.98) — the same pollution-responsive gene cluster identified independently by Rayane
4. **Top Factor 5 genes point to 3 disrupted systems**: immune defence (AZU1), digestion (CTRB1, CPA1), and chitin/exoskeleton metabolism (CHI3L1)
5. **CHI3L1 links to metabolomics pathway findings** — connects to Snehal's amino sugar and nucleotide sugar metabolism pathway (UDP-glucuronic acid → chitin biosynthesis)

---

## How to Reproduce

1. Log into BlueBear HPC
2. `module load bear-apps/2022b/live`
3. `module load MOFA2/1.10.0-foss-2022b-R-4.3.1`
4. Place all input CSV files in working directory
5. Run `Rscript MOFA2_integration.R`
6. All plots and tables will be generated in the same directory
