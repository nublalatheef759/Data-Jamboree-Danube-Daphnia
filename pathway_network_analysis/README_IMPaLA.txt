IMPaLA Multi-Omics Pathway Integration — README

Overview
Integrated pathway over-representation analysis combining transcriptomics (gene list) and metabolomics (KEGG compound IDs) using IMPaLA (http://impala.molgen.mpg.de/). This is the cross-omics validation step — pathways significant here are supported by BOTH gene-level and metabolite-level evidence.

Inputs
- Gene list: meblue_human_symbols.txt (341 human ortholog symbols from Rayane's MEblue WGCNA module)
- Metabolite list: pos_sig_kegg_ids (from Rush, positive mode) + neg_sig_kegg_ids (from Snehal, negative mode), combined into one list
- Gene identifier: HGNC gene symbol
- Metabolite identifier: KEGG compound ID
- Analysis type: Pathway over-representation analysis

Mapping summary
- 145 out of 199 gene identifiers mapped (gene background: 13,434)
- 82 out of 344 metabolite identifiers mapped (metabolite background: 4,359)
- 2,087 pathways tested
- 24 pathways significant at Q_joint < 0.05

Files in this folder
- ORA_results.csv — Full IMPaLA output (2,087 pathways with p-values and q-values for genes, metabolites, and joint test)
- IMPaLA_significant_pathways.xlsx — Curated summary of 24 significant pathways with key genes, metabolite counts, AOP level mapping, and cross-omics validation status

Key findings

Cross-omics validated pathways (genes + metabolites both significant):
1. Metabolism (Reactome) — 50 genes + 16 metabolites, Q = 3.76e-05
2. Transport of small molecules (Reactome) — 23 genes + 6 metabolites, Q = 9.81e-04
3. Metabolism of lipids (Reactome) — 22 genes + 2 metabolites, Q = 8.36e-03
4. Glutathione conjugation (Reactome) — GSTT1, GSTA4 + 1 metabolite, Q = 3.46e-02
5. Aminosugars metabolism (EHMN) — 3 genes + 5 metabolites, Q = 3.46e-02
6. Purine nucleotide degradation (HumanCyc) — NT5E + 4 metabolites, Q = 4.41e-02

Gene-only pathways (biologically important):
- Lysosome (KEGG) — 11 genes including hub genes NPC2, CTSD, GLB1, MAN2B1 (Q = 1.39e-03)
- Angiotensinogen metabolism (Reactome) — hub genes ACE, ANPEP, CTSD, MME (Q = 7.31e-03)
- NRF2 pathway (WikiPathways) — SOD3, GSTA4, GSTT1, master antioxidant regulator (Q = 1.42e-02)
- Neutrophil degranulation (Reactome) — 5 hub genes: ANPEP, CTSD, NPC2, CD36, CAT (Q = 2.76e-02)

Cross-validation with other group members' findings:
- Aminosugars metabolism confirmed — validates Snehal's negative mode metabolomics finding
- Purine nucleotide degradation confirmed — validates Rush's positive mode caffeine/purine finding
- Lysosomal pathway contains 5 hub genes from network analysis — independent confirmation
- NRF2 pathway supports the stress response Key Event (KE1) in the proposed AOP
