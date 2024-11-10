# Thesis Project: Spatial Transcriptomics Analysis
_Short Specialization Degree in Omics Data Analysis, University of Padua_

This repository contains the code and methods used in my thesis project for the short specialization degree in Omics Data Analysis at the University of Padua. The project focuses on spatial transcriptomics, utilizing publicly available datasets from the spatialLIBD package.

##Aim of the project: spatial transcriptomics analysis of Human DLPFC
This project uses spatial transcriptomics data from the Human Dorsolateral Prefrontal Cortex (DLPFC) to analyze gene expression patterns in specific cortical layers and identify spatially variable genes (SVGs). Particularly, the differences between the second layer of the DLPFC and the white matter are analyzed as a proof-of-principle of the methodology. A pdf presentation is also included in this repository, explaining the various steps of the analyses, outputs (tables and plots), and the conclusions.

## Requirements
This analysis relies on the following R packages:

`spatialLIBD` - Data handling and spatial visualization.
`SingleCellExperiment` - Handling single-cell data structures.
`SpatialExperiment` - Defines an S4 class for storing data from spatial -omics experiments.
`nnSVG` - Detection of spatially variable genes.
`PlackettLuce` - Ranking and preference analysis.
`dplyr` and `ggplot2` - Data manipulation and visualization.
`rmarkdown` and `knitr` - Document generation and table rendering.

## Analysis Workflow
1. Data Loading and Preparation
The Visium spatial transcriptomics data for DLPFC samples is loaded and updated using functions from the spatialLIBD package. This step includes converting single-cell data objects into a format suitable for spatial analysis.
2. Data Exploration and Visualization
Spatial Information: Visualizes the distribution of cells in different cortical layers for a specific sample.
Gene Expression: Shows spatial patterns of gene expression for white matter marker genes.
Gene Ranking: Identifies the layer with the highest enrichment of specific genes.
3. Spatial Registration
Single-nucleus RNA-seq data is loaded and used to perform spatial registration, aligning gene expression profiles from single cells to the spatial transcriptomics data. This process computes pseudo-bulking, Bayesian correlation, and t-statistics to enhance spatial resolution.
4. Spatially Variable Gene Analysis
The nnSVG package is used to identify spatially variable genes (SVGs) within specific cortical layers, including white matter and Layer 2. The script demonstrates filtering and downsampling steps to ensure robust SVG analysis and addresses potential issues with duplicate coordinates.
5. Ranking and Comparison of Genes
The PlackettLuce package is used to rank genes based on expression in white matter and Layer 2, generating coefficients for each gene to assess its relative importance. Plots display these rankings for visual comparison between the two layers.

## Example Visualizations
This script includes visualization functions to:
- Display gene expression patterns for specific markers across spatial coordinates.
- Plot layer-specific gene enrichment and SVG patterns for comparison.

## Output
The analysis generates:
- Tables: Summarizing significant genes, spatial registration results, and annotated cell types.
- Plots: Gene expression distribution across layers and ranked gene importance for white matter and Layer 2.
