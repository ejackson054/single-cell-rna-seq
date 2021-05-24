# single-cell-rna-seq
Codes for scRNA-Seq preprocessing and quality control (based on Seurat package, https://satijalab.org/seurat/)

This directory contains a series of scripts for scRNA-Seq data processing
and quality control. They should be run in the following order:

1) Generate_single_cell_matrix.py

Produces counts matrix (gene x cell) for scRNA fastq files.  This script
is a wrapper for kallisto bustools. It runs on one sample at a time; you may
want to write your own wrapper if you have many samples to run.

2) quality_control_pipeline_Part_I.R

Generates set of basic quality control plots in a new folder called 
quality_control. Takes a list of samples as input (each sample must already
have a counts matrix produced by Generate_single_cell_matrix.py).  QC plots
include violinplots broken down by sample (# UMIs, # unique features, %
mitochondrial expression, % ribosomal expression, # doublets), # cells per
sample (table and bar plot), and UMAPs broken down by sample. Uses packages
from Seurat and SCDS (doublet scores).

3) filter_and_integrate_data.R

Filters cells based on user-selected thresholds, performs integration
between different samples, and generates feature plots to assist with
cell type identification.

4) assign_cell_types.R

Allows you to assign cell types based on the clusters and feature plots
created in the previous step.  You manually inspect the clusters and feature
plots, then input a list that contains a cell type for each identified cluster.
The code produces a labeled UMAP plus plots showing the expression of key 
marker genes in each cell type.

5) quality_control_pipeline_Part_II.R

Generate additional quality control plots in existing folder (quality_control). 
QC plots include correlation of gene expression between samples within each 
cell type (as heatmaps and as plots), key QC measures broken down by celltype
(# features, percent mito, percent ribo), expression of sex-linked genes 
(set of Y genes, set of X genes, XIST), UMAP plots for each sample,
and cell type composition for each sample.
