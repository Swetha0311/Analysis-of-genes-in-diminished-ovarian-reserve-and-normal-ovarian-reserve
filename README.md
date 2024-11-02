This repository contains an R script (Deseq2.R) for performing differential expression analysis on RNA-seq count data using the DESeq2 package. This analysis includes data preprocessing, normalization, and visualization of significant genes through heatmaps, volcano plots, and other visual aids.

Requirements

The following R packages are required to run the script:

DESeq2: For differential expression analysis.
ggplot2: For custom plots.
pheatmap: To create heatmaps of gene expression.
EnhancedVolcano: For volcano plot visualization.

Usage

Set Working Directory: Adjust the setwd() function to specify the correct working directory.
Input Data: Place your RNA-seq count data in a .csv file and update the file name in read.table() to your dataset.
The script expects a file named dor_nor.csv, with gene names in the first column and counts in the remaining columns.
Run the Script: Execute the script in R to perform the analysis and generate plots.
Output

The script provides:

Data Summary: Checks for duplicate and missing values.
Differential Expression Analysis: Identification of differentially expressed genes.
Visualizations:
Heatmap of differentially expressed genes.
Volcano plot highlighting significant genes.
