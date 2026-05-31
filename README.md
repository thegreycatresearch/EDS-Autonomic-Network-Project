# 🧬 Autonomic Dysfunction in Ehlers-Danlos Syndrome: A Systems Biology Approach

## 📖 Overview

This project investigates the molecular basis of autonomic dysfunction in Ehlers-Danlos syndrome using integrative transcriptomics, network biology, and pathway enrichment analysis.

We focus on the interaction between:

- Extracellular matrix dysfunction
- Autonomic nervous system signaling
- Ion channel regulation
- Neurovascular pathways

in :contentReference[oaicite:0]{index=0}.

---

## 🎯 Hypothesis

Autonomic dysfunction in EDS arises from a systems-level reprogramming of gene networks linking extracellular matrix abnormalities with neurovascular and autonomic signaling pathways.

---

## 🧪 Methodology

We implemented a reproducible pipeline in R including:

1. GEO dataset retrieval  
2. Quality control and normalization  
3. Differential expression analysis (EDS vs POTS vs controls)  
4. Weighted Gene Co-expression Network Analysis (WGCNA)  
5. Protein interaction network inference  
6. Functional enrichment analysis (GO / KEGG / Reactome)

---

## 🧬 Data Sources

All datasets were obtained from:

:contentReference[oaicite:1]{index=1}

---

## 🔬 Key Analyses

- Differentially expressed genes (DEGs)
- Shared molecular signatures between EDS and POTS
- Gene co-expression modules (WGCNA)
- Hub gene identification
- Neurovascular pathway enrichment

---

## 📊 Outputs

- DEG tables (EDS vs POTS)
- Co-expression modules
- Hub genes
- Enriched biological pathways
- Network models of disease

---

## 🧠 Biological Model

The results support a model where:

- ECM dysfunction → altered mechanotransduction  
- Mechanotransduction → transcriptional reprogramming  
- Reprogramming → autonomic dysregulation  

---

## ⚙️ Reproducibility

This project uses:

- R (≥4.0)
- renv for environment control
- Fully modular scripts

To reproduce:

```r
renv::restore()
source("scripts/01_download_data.R")
source("scripts/02_quality_control.R")
source("scripts/03_differential_expression.R")
source("scripts/04_network_analysis.R")
source("scripts/05_enrichment_integration.R")
