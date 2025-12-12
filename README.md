# Biological Activity of Biomarkers Associated With Metastasis in Osteosarcoma Cell Lines

**08/26/2025**

- Dulce Alejandra Carrillo Carlos [carrillo@lcg.unam.mx](mailto:carrillo@lcg.unam.mx)  

This analysis replicates the following article: [Biological Activity of Biomarkers Associated With Metastasis in Osteosarcoma Cell Lines](https://doi.org/10.1002/cam4.70391).

The data supporting the findings of this study are openly available in **ENA EMBL-EBI** ([https://www.ebi.ac.uk/ena/browser/home](https://www.ebi.ac.uk/ena/browser/home)) under the Project IDs: **PRJEB65286**, **PRJNA698672**, and **PRJNA51801**.

-----

## Problem Statement

Osteosarcoma is the most common cancer among adolescents and young adults, and it represents the most frequent primary malignant bone tumor in this population. It mainly affects long bones during periods of rapid growth, such as the distal femur, proximal tibia, and humerus. Although conventional treatment includes combined surgery and chemotherapy—with significant improvements in survival for patients without disseminated disease—**pulmonary metastasis** remains the leading cause of death. This form of dissemination occurs in approximately 20% of cases at the time of diagnosis and can develop in up to 40% during follow-up. The presence of metastasis is associated with a significantly poorer prognosis and a high recurrence rate, highlighting the urgent need to better understand the molecular mechanisms driving metastatic progression.

At the molecular level, metastasis involves a series of changes in **gene expression**, including pathways related to migration, adhesion, and remodeling of the tumor microenvironment. However, the **molecular biomarkers** associated with these processes are not yet fully characterized, which limits the development of new therapeutic strategies.

Studying gene expression profiles associated with cellular invasion and migration is crucial to:

- Identify **biomarkers** that enable more accurate prognoses.
- Understand the molecular mechanisms of metastasis.

This project aims to **replicate and validate** the findings of the referenced article using a bioinformatics pipeline for RNA-seq analysis and differential expression analysis. Doing so will allow us to explore new hypotheses and generate more robust evidence for future clinical applications.

-----

## Work Plan

| Activity                                                                 | Date             | Responsible     | Deliverable                   |
| ------------------------------------------------------------------------ | ---------------- | --------------- | ----------------------------- |
| Project description                                                      | Aug 28 – Sep 2   | Whole team      | Markdown document             |
| Data download and quality control                                        | Sep 5 – Sep 20   | Andrea           | QC report (FastQC)            |
| Read alignment and quantification                                        | Sep 21 – Sep 30  | Antonio          | Files, count table            |
| Differential expression analysis (DESeq2)                                | Oct 1 – Oct 15   | Alejandra        | Script, DEG tables            |
| Functional annotation                                                    | Oct 16 – Oct 25  | Alejandra        | Markdown report               |
| Results integration and visualization                                    | Oct 26 – Nov 5   | Antonio          | Figures in repository         |
| Writing of results report                                                | Nov 6 – Nov 15   | Andrea           | Markdown document             |
| Project presentation                                                     | December         | Whole team      | GitHub repository (release)   |

-----

## Methodology

### Research Question

1. Which genes are differentially expressed between normal osteoblastic cells and metastatic osteosarcoma cell lines?

The conditions of the samples are the following, OS means Osteosarcoma and CTL means control:

| **Run Accession** | **Sample Accession** | **Sample Title** |
|:-----------------:|:--------------------:|:----------------:|
| `ERR13451008` | `SAMEA114285858` | OS |
| `ERR13448321` | `SAMEA114285854` | CTL |
| `ERR13451021` | `SAMEA114285857` | OS |
| `ERR13448329` | `SAMEA114285856` | OS |
| `ERR13451054` | `SAMEA114285859` | OS |
| `ERR13447502` | `SAMEA114285855` | CTL |

### Proposed Pipeline

1. **Data retrieval and download** (ENA; projects PRJEB65286, PRJNA698672, PRJNA51801)  
2. **FASTQ quality control** (FastQC, Trim Galore)  
3. **Read alignment** (STAR)  
4. **Count matrix generation**  
5. **Differential expression analysis** (DESeq2)  
6. **Functional annotation**  
7. **Results visualization**  
8. **Integration and discussion of results**  

-----

## Expected Results

- Identification of differentially expressed genes between normal cells and metastatic osteosarcoma.  
- Detection of relevant biomarkers involved in metastasis.  
- Reproduction of key results reported in the original article, including figures.  
- Creation of a reproducible GitHub repository.
