Network Analysis — README

Overview
Protein-protein interaction (PPI) network analysis and hub gene identification from the MEblue WGCNA module (741 genes, major pollution response cluster). Hub genes were used to construct a proposed Adverse Outcome Pathway (AOP) for Danube pollution effects on Daphnia magna. Hub genes were also checked against known human disease associations using GeneCards and OMIM.

Input
- meblue_human_symbols.txt — 341 human ortholog gene symbols from Rayane's MEblue WGCNA module (mapped from Daphnia magna via D. pulex intermediary). This was the input for STRING.

Tools used
1. STRING (https://string-db.org) — PPI network construction. Organism: Homo sapiens (human orthologs). Full STRING network, medium confidence (0.4). Result: 187 nodes, 240 edges.
2. Cytoscape — Network visualisation and analysis. Imported STRING TSV export.
3. cytoHubba (Cytoscape plugin) — Hub gene identification using three graph-theoretic centrality methods: MCC (Maximal Clique Centrality), Degree, and Betweenness.
4. DiagrammeR (R package) — AOP diagram construction.
5. GeneCards (https://www.genecards.org) and OMIM (https://omim.org) — Human disease association lookup for hub genes.

Files in this folder

| File | Description |
|------|-------------|
| string_vector_graphic.svg | Full STRING PPI network (187 nodes, 240 edges) |
| string_interactions_short_tsv_MCC_top10.svg | cytoHubba top 10 hub genes ranked by MCC |
| string_interactions_short_tsv_Degree_top10.svg | cytoHubba top 10 hub genes ranked by Degree |
| string_interactions_short_tsv_Betweenness_top10.svg | cytoHubba top 10 hub genes ranked by Betweenness |
| hub_gene_comparison_table.xlsx | Comparison table of hub genes across all 3 methods, with biological roles and AOP level mapping |
| hub_gene_disease_associations.xlsx | Human disease associations for each hub gene, sourced from GeneCards and OMIM |
| AOP_diagram.svg | Proposed Adverse Outcome Pathway diagram with hub genes at each Key Event level |

Key findings

Top hub genes (appearing in all 3 centrality methods)
- ITGB1 — integrin beta-1, cell adhesion and signalling (MIE)
- ANPEP — aminopeptidase, protein processing (KE2)
- CD36 — lipid/fatty acid uptake receptor (MIE)
- ACE — angiotensin-converting enzyme, also a universal pollution responder from WGCNA (KE1)

Strong hub genes (appearing in 2/3 methods)
- ABCA1 — cholesterol efflux transporter (KE1)
- ABCG2 — xenobiotic efflux pump (KE1)
- HSP90AA1 — heat shock protein, stress response (KE1)
- CTSD — lysosomal protease (KE2)
- CD81 — tetraspanin, membrane signalling (MIE)

Human disease relevance
Hub genes map to known human diseases, supporting conservation of these pathways:
- ACE — hypertension, cardiovascular disease, diabetic nephropathy (target of ACE inhibitor drugs)
- CD36 — coronary heart disease susceptibility (OMIM)
- ABCA1 — Tangier disease, coronary artery disease
- ABCG2 — gout, hyperuricemia (uric acid transport)
- CTSD — neuronal ceroid lipofuscinosis (Batten disease, lysosomal neurodegeneration)
- CD81 — immunodeficiency (CVID6), hepatitis C susceptibility
- HSP90AA1 — glioma susceptibility, hematologic cancer

Proposed AOP structure
Stressor (Danube pollution) → MIE: membrane receptor disruption (ITGB1, CD36) → KE1: stress response and xenobiotic efflux (HSP90AA1, ACE, ABCG2) → KE2: lipid and lysosomal disruption (NPC2, ANPEP, CTSD) → KE3: suppressed growth and reproduction (~10:1 down:up ratio) → AO: reduced fitness and survival of Daphnia magna

Connection to other analyses
- ACE was independently identified as a universal pollution responder in Rayane's WGCNA (DEG in >=15/24 contrasts), confirming its importance through graph theory
- Hub genes are dominated by membrane/lipid/lysosomal biology, consistent with Rayane's finding that lipid metabolism (FABP12, NPC2, LIPC) dominates the MEblue module
- Diuron has the strongest single-chemical correlation with MEblue (r = +0.585)
- D11 (1,570 DEGs) and D12 (1,605 DEGs) show the most severe transcriptional disruption
- IMPaLA cross-omics integration confirmed aminosugars metabolism (Snehal's finding) and purine nucleotide degradation (Rush's finding) are significant when combining genes and metabolites
