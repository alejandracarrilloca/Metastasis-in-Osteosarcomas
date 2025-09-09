# Biological Activity of Biomarkers Associated With Metastasis in Osteosarcoma Cell Lines


## Reference Organism

Osteosarcoma is a malignant bone tumor that affects humans. In this project, the RNA-seq data come from human metastatic osteosarcoma cell lines and normal human osteoblastic cells.

The analysis was performed using the human genome (Homo sapiens), GRCh38.p12 (RefSeq: GCF\_000001405.33, available at NCBI) as the reference. This genome is currently the standard for gene expression studies, as it contains the most up-to-date annotation of genes, isoforms, and regulatory regions.

The RNA-seq reads were aligned and quantified against this genome using Kallisto (Pachter Lab), which enables fast and accurate quantification of transcripts. The use of GRCh38.p12 ensures comparability of results with other human cancer studies and allows mapping of genes involved in invasion, migration, and metastasis processes.

## Objectives

- Analyze gene expression profiles in normal osteoblastic cells and metastatic osteosarcoma cell lines using RNA-seq data.
- Identify differentially expressed genes that establish clear differences between normal and tumor conditions, with emphasis on those involved in invasion and metastasis processes.
- Detect relevant biomarkers associated with metastatic progression, which may have potential value as prognostic predictors or therapeutic candidates.

## Expected Results

Differential expression analysis is expected to reveal a robust set of genes characteristic of osteosarcoma with metastatic potential, enriched in tumor progression pathways, including ECM–receptor interaction, focal adhesion, MAPK signaling, and actin cytoskeleton regulation.

From these genes, a functional classification will be performed, highlighting biological processes directly linked to metastasis, such as cell adhesion, migration, extracellular matrix remodeling, and tumor signaling. In protein–protein interaction networks, hubs are expected to emerge, such as collagens and fibronectin (matrix components), integrins (adhesion receptors), growth factors (FGF2), and proteases (MMP-2), consistent with their previously described role in other cancer types.

Following this analysis, the identification of potential biomarkers with clinical utility for prognosis is anticipated, and eventually, these could guide targeted therapeutic strategies in metastatic osteosarcoma.

Additionally, partial validation and reproduction of findings reported in the literature are expected. Finally, the creation of a GitHub repository is planned to document the scripts, tables, and key figures of the analysis, facilitating reproducibility of the study.

## BioProject ID

The data supporting the findings of this study are openly available in ENA EMBL-EBI under the Project IDs: PRJEB65286, PRJNA698672, and PRJNA518013.


