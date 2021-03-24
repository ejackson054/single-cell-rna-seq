#This script will:  1) Filter cells based on criteria from Quality Control Pipeline #1, 2) Integrate data between samples, 
#3) Add in your metadata to the Seurat object (minimum batch and karyotype), 4) Make a UMAP, 5) Generate feature plots to help assign cell types

library(BUSpaRse)
library(ggplot2)
library(magrittr)
library(data.table)
library(Seurat)
library(DropletUtils)
library(Matrix)
theme_set(theme_bw())
library(reticulate)
library(dplyr)
library(scds)
library(scater)
library(rsvd)
library(Rtsne)
library(cowplot)
library(hash)
library(reshape)

setwd("/lab/page_human_data/emily/single_cell/")

###SET QC THRESHOLDS###

min_features <- 1000
max_mito <- 15
max_doublet <- 0.8

###CHOOSE SAMPLES###   these should be the names of samples you want to include from the bus_output folder (excluding any outliers!)
samples_list <- c("L16_2576_S1","L16_2577_S2")

###READ IN GENE METADATA###
meta <- read.csv("/lab/Page_lab-users/Alex/gtex/index/gencode.v24.annotation.basic_ccds_nopar.gene_tx_annotable.txt",sep="\t")
meta <- distinct(meta, gene_name, .keep_all = TRUE)
rownames(meta) <- meta$gene_name

genes_of_interest <- c()

for (i in 1:length(rownames(meta))) {
  gene <- rownames(meta)[i]
  if (meta[gene,"gene_type"]=="protein_coding") {
    genes_of_interest <- c(genes_of_interest,gene)
  }
  else if (gene %in% c("XIST","TXLNGY", "PRKY")) {
    genes_of_interest <- c(genes_of_interest,gene)
  }
}


###GENERATE FILTERED SEURAT OBJECT FOR EACH SAMPLE###
object_list <- c()

for (sample in samples_list) {
  
  #Read in file
  directory <- sprintf("bus_output/%s/genecount", sample)
  res_mat <- read_count_output(directory, name = "gene", tcc = FALSE)
  
  #Filter to protein-coding genes only
  res_mat <- res_mat[genes_of_interest,]
  
  #Make object
  object <- CreateSeuratObject(counts=res_mat, min.cells = 3, min.features=min_features, project=sample) 
  
  #Find fractions mitochondrial and platelet
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  object[["percent.platelet"]] <- PercentageFeatureSet(object, features = c("PPBP","PF4"),assay="RNA") #"GP1BB","CAVIN2","GNG11"
  object[["percent.ribo"]] <- PercentageFeatureSet(object, pattern = "^RP[LS]")
  
  #Add doublet scores
  sce <- as.SingleCellExperiment(object)
  sce <- cxds(sce,retRes = TRUE)  #coexpression
  sce <- bcds(sce,retRes = TRUE,verb=TRUE) #artificial doublets
  sce <- cxds_bcds_hybrid(sce)
  
  object@meta.data$hybrid_score <- sce$hybrid_score
  object@meta.data$bcds_score <- sce$bcds_score
  object@meta.data$cxds_score <- sce$cxds_score

  #Filter and normalize object
  object <- subset(object, subset = percent.mt < max_mito & bcds_score <max_doublet & percent.platelet < 0.01)
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object<- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
  
  #Add to list
  object_list <- c(object_list,object)
  
}

###INTEGRATE DATA###
seu.anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:30, reference = c(1)) #takes around 15 min for 4 samples
seu <- IntegrateData(anchorset = seu.anchors, dims = 1:30)
DefaultAssay(seu) <- "integrated"

###READ IN SAMPLE METADATA### 
meta <- read.csv(file="metadata.csv")

hash.batch <- hash(meta$Sample_ID, meta$Batch)
hash.kar <- hash(meta$Sample_ID, meta$Karyotype)

batches <- c()
karyotypes <- c()

for (i in rownames(seu@meta.data)) {
  batch <- hash.batch[[seu@meta.data[i,"orig.ident"]]]
  batches <- c(batches, as.character(batch))
  karyotype <- hash.kar[[seu@meta.data[i,"orig.ident"]]]
  karyotypes <- c(karyotypes, as.character(karyotype))
}

seu <- AddMetaData(seu, batches, col.name = "Batch")
seu <- AddMetaData(seu, karyotypes, col.name = "Karyotype")


###ADD CELL CYCLE INFO###    https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette.html
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seu$CC.Difference <- seu$S.Score - seu$G2M.Score

###SAVE OBJECT###
if (!dir.exists("integrated_results")) {
  system(sprintf("mkdir integrated_results"))
}

saveRDS(seu, file = "integrated_results/integrated_data.rds")

###SCALE DATA AND RUN PCA###
set.seed(8)

if (length(unique(seu@meta.data$Batch))>1) {
  seu <- ScaleData(seu,verbose=FALSE,vars.to.regress = c("percent.mt", "Batch", "CC.Difference"))
} else {
  seu <- ScaleData(seu,verbose=FALSE,vars.to.regress = c("percent.mt", "CC.Difference"))
}

  
set.seed(8)
seu <- RunPCA(seu, verbose = FALSE, npcs = 100)

#Find clusters and run UMAP
set.seed(8)
seu <- FindNeighbors(seu, reduction = "pca", dims=1:75,  nn.eps = 0.5) #lower = more exact, slower:  default = 0
set.seed(8)
seu <- FindClusters(seu, resolution = 1.5, graph.name = 'integrated_snn')
set.seed(8)
seu <- RunUMAP(seu, min.dist = 0.3, dims=1:75, metric="euclidean", n.neighbors=20L)  


###MAKE INTEGRATED UMAP###
png("integrated_results/umap_integrated.png", width = 800, height = 600, units = "px")
DimPlot(seu, reduction = "umap", label=TRUE, repel=TRUE)
dev.off()

###MAKE FEATURE PLOTS FOR CELL TYPE ASSIGNMENT###
DefaultAssay(seu) <- "RNA"

if (!dir.exists("integrated_results/feature_plots")) {
  system(sprintf("mkdir integrated_results/feature_plots"))
}

png("integrated_results/feature_plots/B_Monocytes.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("MS4A1", "CD14", "FCGR3A"))
dev.off()

png("integrated_results/feature_plots/Dendritic_MAIT.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("FCER1A","UGCG","SLC4A10"))
dev.off()

png("integrated_results/feature_plots/NK.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("CD3D","GNLY","SELL"))
dev.off()

png("integrated_results/feature_plots/CD8_T.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("CD8A","CCR7","GZMB"))
dev.off()

png("integrated_results/feature_plots/CD4_T.png", width = 1500, height = 500, units = "px")
FeaturePlot(seu, ncol = 3, features = c("CD4","CCR7","FOXP3"))
dev.off()


###SAVE FINAL FILE###
saveRDS(seu, file = "integrated_results/integrated_data_processed.rds")

