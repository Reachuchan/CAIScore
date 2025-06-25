# 1. Introduction
   **CAIScore (Cellular Activity Integration Score)** is a **generalizable and biologically informed computational framework** designed to **quantify the activity of any user-defined gene set** across **bulk, single-cell, and spatial transcriptomic datasets**. Originally developed to assess DNA damage repair (DDR) dynamics in hepatocellular carcinoma (HCC), **CAIScore integrates gene set enrichment with transcriptomic topology** to provide a **unified, scalable metric** that reflects **both gene set activity and its microenvironmental impact at high resolution**.
   **CAIScore is built on a modular scoring architecture** that allows researchers to input **any curated gene set**—ranging from hallmark oncogenic pathways to immune signatures, metabolic programs, or tissue-specific regulators—thereby enabling **flexible application across a wide array of biological questions and disease contexts**. In our study, **CAIScore effectively captured DDR-driven heterogeneity and immunosuppressive remodeling in HCC**, **outperforming conventional scoring methods** in **prognostic prediction and immune landscape stratification**. Notably, **CAIScore maintained robust performance across multiple independent datasets and cancer types**, underscoring its **broad applicability and translational potential**.
   Beyond static scoring, **CAIScore offers mechanistic insight into the dynamic tumor ecosystem** by integrating **gene activity with cellular composition, intercellular communication, and spatial architecture**. The algorithm has been **systematically validated in multi-omics settings**—including bulk RNA-seq, single-cell RNA-seq, and spatial transcriptomics—highlighting its **capacity to dissect functional states, lineage plasticity, and immune evasion mechanisms at scale**.
   This repository provides an **open-source implementation of CAIScore**, complete with **customizable input options, example datasets, and reproducible workflows**. By enabling **biologically meaningful and context-aware activity quantification of any gene set**, **CAIScore serves as a powerful tool** for **uncovering key regulatory programs and therapeutic vulnerabilities across diverse biological systems**.
## 1.1 Overview of the CAIScore Algorithm and Its Multi-Omics Applications in DDR Activation
![CAIScore](https://github.com/Reachuchan/Figure/blob/main/GitHub%20Figure.png)
## 2. Installation
```R
if (!require("devtools")) 
  install.packages("devtools")
devtools::install_github("Reachuchan/CAIScore")
```
## 3. Using the Code
```R
library(UCell)
library(AUCell)
library(singscore)
library(GSVA)
library(data.table)
library(tibble)
library(CAIScore)

## Bulk RNA-seq
DDR_Genes = DDR_Genes

expr = bulk_exp

score_result <- CAIScore(
  expr = expr,              # Expression Matrix (Bulk-RNA)
  geneset = DDR_Genes,
  scaling = "minmax",       # Options: "zscore", "minmax", "none"
  summary_method = "sum",   # Options: "mean", "sum", "median"
  auc_nCores = 1,           # Number of threads for AUCell
  kcdf_type = "Gaussian",   # For normalized data: use "Gaussian";For count data: use "Poisson"
  #AMS_nbin = 10             #(default:24);#Useful to reduce this (e.g. to 10) when sample size is small to avoid cut_number errors.
)

## single-cell RNA-seq/spatial transcriptome RNA-seq

DDR_Genes = DDR_Genes

expr = scRNA #stRNA

score_result <- CAIScore(
  expr = expr,              #seurat scRNA/stRNA
  geneset = DDR_Genes,
  scaling = "minmax",       # Options: "zscore", "minmax", "none"
  summary_method = "sum",   # Options: "mean", "sum", "median"
  auc_nCores = 1,           # Number of threads for AUCell
  kcdf_type = "Gaussian",   # For normalized data: use "Gaussian";For count data: use "Poisson"
  AMS_nbin = 10             #(default:24);#Useful to reduce this (e.g. to 10) when sample size is small to avoid cut_number errors.
)
```
