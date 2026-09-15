# Comparative Analysis of scRNA-Seq Data under Drought Stress Conditions
This repository contains reproducible workflows and pipelines for the quality control, integration, downstream analysis, and network construction of rice scRNA-seq datasets. The workflow demonstrates computational approaches to dissect cellular heterogeneity and gene regulatory mechanisms underlying the response of rice to drought stress.

## Key Capabilities & Workflows
* **Data Quality Control & Preprocessing:** Adaptive filtering of low-quality cells based on UMI counts, feature counts, and organellar (mitochondrial/chloroplast) gene expression thresholds.
* **Dataset Integration:** Implementation of integration algorithms (`McInnes et al., 2018`) to harmonize the normal and drought datasets by reducing unwanted sources of variation and facilitating the alignment of shared cell populations across conditions, while preserving biologically meaningful, cell type-specific transcriptional differences associated with drought stress.
* **Cell-Type Annotation:** Cell types were annotated using MICI (`Wang et al., 2021`), followed by validation based on the expression patterns of literature-curated marker genes.
* **Condition-Specific DEG Analysis (Normal vs. Drought):** Differential expression testing between experimental conditions within specific cell types using `FindMarkers` to identify stress-responsive genes and cell-type-specific transcriptional shifts.
* **Network Modeling:**
  * **hdWGCNA:** High-dimensional Weighted Gene Co-expression Network Analysis to discover co-expressed gene modules specific to cell types (`Morabito et al., 2023`).
  * **SCENIC:** Gene Regulatory Network (GRN) inference and regulon activity scoring (`Aibar et al., 2017`) were used to identify cell-type-specific transcriptional regulators and their target genes associated with the rice response to drought stress. My pipeline utilizes **GRNBoost2** for fast, scalable network inference in the initial stage, followed by downstream regulon analysis in R.

<img width="800" height="450" alt="Figure 1  workflow_scRNA_PP" src="https://github.com/user-attachments/assets/40d6c309-aae6-4c76-bfaa-1a15bd8c7f3f" />

## Software and Computational Tools
* **Languages:** R, Python
* **Core Libraries:** Seurat, Bioconductor, Tidyverse, hdWGCNA, SCENIC

## Repository Structure
```text
├── data/               
├── src/
│   ├── 01_QC_and_Normalization.R    # Quality control and SCTransform workflows
│   ├── 02_Data_Integration.R        # Data integration using Harmony
│   ├── 03_Clustering_Annotation.R   # Cell typing
│   ├── 04_Condition_DEG_Analysis.R  # Differential expression between normal vs. drought
│   ├── 05_hdWGCNA_Analysis.R        # Co-expression network modeling
│   └── 06_SCENIC_GRN_Inference.R    # Regulon activity analysis
└── README.md
```

## Acknowledgments
Parts of the methodology implemented in this repository were utilized in research presented at:
* **International Competition on Science, Technology, Engineering, and Mathematics (ICSTEM 2025), Malaysia** (Awarded Double Gold Medal and Special Award from the Indonesian Invention and Innovation Promotion Association(INNOPA)).
* **The 6th National Rice Research Conference (2025), Thailand** (Oral Presentation: Evaluation of Marker-Based Cell Identity Index in plant single-cell analysis).
* **The Graduate Innovation and Exhibition: National and International (GIENI 2026), Thailand** (Awarded Gold Medal).

---
**Author:** Saranyu Maensatit  
**Contact:** bksaranyu@gmail.com
