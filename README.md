# Biological Activity of Biomarkers Associated With Metastasis in Osteosarcoma Cell Lines

**08/26/2025**

- Dulce Alejandra Carrillo Carlos [carrillo@lcg.unam.mx](mailto:carrillo@lcg.unam.mx)  
- Addiel Antonio Platas Renteral [addielpr@lcg.unam.mx](mailto:addielpr@lcg.unam.mx)  
- Andrea Villarruel García [andreavg@lcg.unam.mx](mailto:andreavg@lcg.unam.mx)

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

## Results Visualization

This repository includes comprehensive visualization scripts to analyze and present the differential expression results:

### Available Visualizations

1. **Volcano Plot** - Shows differential expression significance vs fold change
2. **PCA Plot** - Displays sample clustering and variance decomposition  
3. **Expression Heatmaps** - Visualizes expression patterns of top differentially expressed genes
4. **Pathway Diagrams** - Illustrates enriched biological pathways and metastasis processes

### Quick Start

To generate all visualizations:

```bash
# Install dependencies
pip install -r requirements.txt

# Generate all plots using sample data
python scripts/integrate_results.py

# Or run individual scripts
python scripts/volcano_plot.py -i data/deseq2_results.csv -o figures
python scripts/pca_plot.py -i data/pca_results.csv -o figures  
python scripts/heatmap.py -e data/expression_matrix.csv -d data/deseq2_results.csv -o figures
python scripts/pathway_diagrams.py -d data/deseq2_results.csv -o figures
```

### For Real Analysis Data

When you have your own DESeq2 results, simply replace the sample data files:

```bash
# Replace sample data with your results
cp your_deseq2_results.csv data/deseq2_results.csv
cp your_expression_matrix.csv data/expression_matrix.csv

# Run the complete visualization pipeline
python scripts/integrate_results.py
```

### Output Files

The visualization scripts generate both PNG and PDF formats:
- `volcano_plot.png/pdf` - Differential expression volcano plot
- `pca_plot.png/pdf` - PCA analysis with variance explanation
- `heatmap.png/pdf` - Expression heatmap of top DE genes
- `clustered_heatmap.png/pdf` - Hierarchically clustered heatmap
- `pathway_enrichment.png/pdf` - Pathway enrichment analysis
- `metastasis_pathway_diagram.png/pdf` - Conceptual metastasis pathway

All figures are saved to the `figures/` directory by default.
