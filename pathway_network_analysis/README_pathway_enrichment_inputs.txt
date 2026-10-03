PATHWAY ENRICHMENT INPUTS
=========================

Use these files for GO/KEGG/Reactome enrichment analysis.

KEY FILES:
- A_pooled_10x_vs_Control.csv : full ranked DEG list, all 17,881 genes
                                (use for GSEA-style methods needing ranked lists)
- B_pooled_1x_vs_Control.csv  : same for low-dose contrast
- universal_responders_anchored.csv : 29 universal pollution responders
                                       (DEG in 15+ of 24 conditions)
- top50_pooled_10x_DEGs_annotated.csv : top 50 DEGs with human ortholog mapped
- meblue_genes_annotated.csv : 741 genes in WGCNA MEblue module
                                (Diuron-linked detox cluster)

ORTHOLOG MAPPING:
- ortholog_mapping_anchored.csv : full Dma -> dpx -> dme -> hsa table
  Anchored on D. pulex (closest species), then propagated.
  Resolves 1:many by primary (first-listed) ortholog.

DEG THRESHOLD USED: padj < 0.05 AND |log2FC| > 1

