# Data Directory

This folder is intentionally empty. The original data files are **not included** in this repository due to intellectual property restrictions.

## Required Data Files

To reproduce the analyses, you need the following files (provided by the module organisers):

| File | Description |
|------|-------------|
| `sample_sheet.csv` | 150 samples: 6 Control + 72 × 1x + 72 × 10x across 12 Danube sites |
| `rna_raw_counts.csv` | Raw RNA-seq count matrix (17,881 genes × 150 samples) |
| `rna_norm_counts.csv` | DESeq2-normalised counts |
| `rna_vst_counts.csv` | Variance-stabilising-transformed counts |
| `polar_pos_pqn_imputed_glog.csv` | Positive-mode metabolomics (PQN + KNN + glog) |
| `polar_neg_pqn_imputed_glog.csv` | Negative-mode metabolomics (PQN + KNN + glog) |
| `water_chemicals.tsv` | 91 water chemicals measured at 12 Danube sites |
| `rna_dma_to_hsa_mappings.tsv` | Ortholog mapping: *D. magna* → Human |
| `rna_dma_to_dme_mappings.tsv` | Ortholog mapping: *D. magna* → *Drosophila* |
| `rna_dma_to_dpx_mappings.tsv` | Ortholog mapping: *D. magna* → *D. pulex* |
| `polar_pos_pkl_to_kegg_annotations.tsv` | Positive-mode KEGG annotations |
| `polar_pos_pkl_to_hmdb_annotations.tsv` | Positive-mode HMDB annotations |
| `polar_neg_pkl_to_kegg_annotations.tsv` | Negative-mode KEGG annotations |
| `polar_neg_pkl_to_hmdb_annotations.tsv` | Negative-mode HMDB annotations |

Place all data files in this directory before running the analysis scripts.
