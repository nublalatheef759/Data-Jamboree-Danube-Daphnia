# Multi-Omics Ecotoxicology of *Daphnia magna* Along the Danube River

[![R](https://img.shields.io/badge/R-≥4.3-blue)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-DESeq2%20|%20MOFA2-green)](https://www.bioconductor.org/)
[![License](https://img.shields.io/badge/License-Educational-yellow)](LICENSE)
[![Multi-Omics](https://img.shields.io/badge/Multi--Omics-Transcriptomics%20|%20Metabolomics-orange)](.)

**Data Jamboree Group Project — MSc Bioinformatics, University of Birmingham (2026)**

## Research Question

> How do different chemical pollution mixtures along the Danube River affect the biology of *Daphnia magna*, and which pollutants and biological pathways have the greatest impact?

## Project Overview

This project integrates transcriptomics and metabolomics data from *Daphnia magna* exposed to water samples collected at 12 sites along the Danube River at two concentration levels (1x and 10x). The study applies a multi-omics approach to identify pollution-responsive genes, metabolic pathways, and the chemical drivers behind the observed biological effects.

### Experimental Design

- **Organism**: *Daphnia magna* (freshwater crustacean model)
- **Sites**: 12 locations along the Danube (D01–D12)
- **Doses**: Control, 1x (environmental), 10x (concentrated)
- **Samples**: 150 total (6 control + 72 per dose, biological replicates per site)
- **Omics layers**: RNA-seq (17,881 genes), polar metabolomics positive mode, polar metabolomics negative mode
- **Chemical data**: 91 water chemicals measured at each site

## Key Findings

1. **Two distinct pollution regimes** along the Danube:
   - **Upstream–middle (D02–D04, D09)**: agricultural and personal-care chemicals (Diuron, DEET, Terbuthylazine)
   - **Downstream (D11–D12)**: municipal wastewater signature (Metformin, artificial sweeteners, pharmaceuticals)

2. **Wastewater sites D11/D12 show the largest biological response** — 1,570 and 1,605 DEGs respectively at 1x concentration

3. **D12 carries a unique pesticide pulse** (Chlorpropham, Carbaryl, Diazinon, Simazine, Chlorophene) found at no other site

4. **29 universal pollution responders** appear across 15+ of 24 site–dose conditions, clustering into:
   - Lipid metabolism (NPC2, FABP12, LIPC)
   - Phase II detoxification (SULT1A5)
   - Digestive proteases (suppressed)

5. **MOFA2 Factor 5** is the shared multi-omics signal (7.5% RNA + 4.6% neg + 5.1% pos variance), correlating with WGCNA MEblue module (r = 0.98)

6. **72 significant hub-gene × chemical correlations** identified, strongest: CD36 × Mecoprop (r = 0.830), ABCG2 × 2,4-DNP (r = 0.832)

7. **24 significant multi-omics pathways** via IMPaLA (Q < 0.05), with cross-omics validation of aminosugar and purine degradation pathways

## Repository Structure

```
Data-Jamboree-Danube-Daphnia/
├── README.md                              # This file
├── .gitignore
├── data/                                  # Empty — see data/README.md for required files
│   └── README.md
│
├── transcriptomics/
│   └── transcriptomics.R                 #   DESeq2, PCA, QC, WGCNA, volcano plots, ortholog mapping
│
├── metabolomics/                          # Positive- and negative-mode metabolomics
│   ├── metabolomics_positive_EDA.R       #   Pos mode: QC, EDA, ANOVA, PERMANOVA, volcano plots
│   ├── metabolomics_positive_integration.R #  Pos mode: pathway annotation & team handoff
│   ├── metabolomics_negative_EDA.R       #   Neg mode: PCA, ANOVA, PERMANOVA, heatmap, annotation
│   └── metabolomics_negative_integration.R #  Neg mode: integration deliverables (169 sig features)
│
├── integration/                           # Multi-omics integration
│   ├── MOFA2_integration.R               #   MOFA2: 10 factors, 3 views, 142 shared samples
│   └── README.md                         #   Variance explained & key results
│
├── pathway_network_analysis/              # Pathway enrichment & network analysis
│   ├── pathway_network_analysis.R        #   DESeq2 per-site, orthologs, STRING, AOP, correlations
│   ├── README_network_analysis.txt       #   STRING/Cytoscape/cytoHubba hub gene results
│   ├── README_chemical_correlation.txt   #   72 significant gene × chemical correlations
│   └── README_IMPaLA.txt                 #   24 significant multi-omics pathways
│
├── README_transcriptomics_results.txt     # Shared DEG counts, WGCNA modules, D11/D12 chemistry
├── README_integration_inputs.txt          # Module eigengenes, gene assignments, correlations
├── README_pathway_enrichment_inputs.txt   # Pooled DEG lists, universal responders, MEblue genes
└── README_presentation_figures.txt        # Guide to key presentation figures
```

## Methods Summary

| Analysis | Tool / Method | Script |
|----------|--------------|--------|
| Differential gene expression | DESeq2 (padj < 0.05, \|log2FC\| > 1) | `transcriptomics/transcriptomics.R` |
| Gene co-expression | WGCNA (17 modules, top 5,000 genes) | `transcriptomics/transcriptomics.R` |
| Metabolomics QC & normalisation | PQN + KNN imputation + glog transform | `metabolomics/metabolomics_positive_EDA.R` |
| Metabolomics differential | ANOVA + PERMANOVA + volcano plots | `metabolomics/` (both modes) |
| Multi-omics integration | MOFA2 (10 factors, 3 views) | `integration/MOFA2_integration.R` |
| Ortholog mapping | OrthoDB via D. pulex → Drosophila → Human | `pathway_network_analysis/pathway_network_analysis.R` |
| Pathway enrichment | WebGestalt (ORA), iDEP, IMPaLA | `pathway_network_analysis/pathway_network_analysis.R` |
| Network analysis | STRING + Cytoscape + cytoHubba (MCC, Degree, Betweenness) | `pathway_network_analysis/pathway_network_analysis.R` |
| Chemical–gene correlation | Spearman correlation (padj < 0.05) | `pathway_network_analysis/pathway_network_analysis.R` |

## Data Availability

The original data files are **not included** in this repository due to intellectual property restrictions from the module organisers. See `data/README.md` for a list of required input files.

## Software Requirements

- R (>= 4.3)
- Key R packages: DESeq2, WGCNA, pheatmap, ggplot2, corrplot, matrixStats, patchwork
- MOFA2 (via Bioconductor)
- External tools: WebGestalt, iDEP, IMPaLA, STRING, Cytoscape + cytoHubba, MetaboAnalyst
