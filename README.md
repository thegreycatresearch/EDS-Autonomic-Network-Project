# EDS + Dysautonomia Systems Biology

## Overview

This repository contains a systems biology pipeline analyzing transcriptomic datasets associated with **Ehlers-Danlos Syndrome (EDS)** and **dysautonomia (POTS-like phenotypes)**.

The goal is to identify shared molecular mechanisms linking:

- Extracellular matrix (ECM) dysfunction  
- Vascular instability  
- Autonomic dysregulation  

and to build an integrated mechanistic model.

---

## Study Design

The pipeline includes:

1. GEO dataset acquisition  
2. Quality control & normalization  
3. Differential gene expression analysis  
4. Co-expression network analysis (WGCNA-inspired)  
5. Functional enrichment (GO analysis)  
6. Integrated systems-level modeling  
7. Publication-ready figure generation  

---

## Key Figures

- **Figure 1:** Analysis workflow  
- **Figure 2:** Differential gene signatures (EDS vs POTS)  
- **Figure 3:** Co-expression network structure  
- **Figure 4:** Functional enrichment analysis  
- **Figure 5:** Integrated mechanistic model  
- **Graphical Abstract:** System-level disease model  

---

## Biological Hypothesis

We propose a unified model where:

> ECM dysfunction → vascular instability → autonomic dysregulation

forms a continuous biological axis underlying dysautonomia in connective tissue disorders.

---

## Requirements

```r
R version >= 4.2

Packages:
- GEOquery
- Biobase
- ggplot2
- dplyr
- igraph
- patchwork
- pheatmap
```
---

# How to run

## Step 1: Download data
source("scripts/01_download_data.R")

## Step 2: QC
source("scripts/02_quality_control.R")

## Step 3: Analysis pipeline
source("scripts/03_differential_expression.R")

## Step 4: Figures
source("scripts/06_generate_figures.R")

---

## Output

All outputs are saved in:

```r
results/
```
---

## Citation

If you use this pipeline, please cite:

|"Systems biology analysis of dysautonomia in connective tissue disorders"

---

## Author

Systems Biology Project — EDS / Dysautonomia Transcriptomics
Bianca Stazzone
