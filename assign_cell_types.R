#Script takes a list of assigned cell types, and creates a) new labeled UMAP, b) violin plots showing the expression of 
#key marker genes in each cell type, c) dot plot showing the expression of key marker genes in each cell type

#In simplest form, script can be run after defining marker genes (line 24) and assigned cell types for each cluster (line 31).
#If desired, colors on UMAP and order of cell types in dot plot can later be tweaked using the code in lines 79 and 100, respectively.


library(ggplot2)
library(magrittr)
library(data.table)
library(Seurat)
library(Matrix)
library(dplyr)
library(hash)
library(cowplot)
library(scales)


###SET VARIABLES###
base_directory <- "/lab/page_human_data/emily/single_cell/integrated_results"  #folder that contains results from filter_and_integrate_data.R
rds_file <- "integrated_data_processed.rds"   #RDS file created by filter_and_integrate_data.R


###SET MARKER GENES###        
#List of genes used in feature plots from filter_and_integrate_data.R
marker_genes <- c("MS4A1","CD14","FCGR3A", "FCER1A", "UGCG",  
                  "GNLY",  "CD3D", "FOXP3", "SLC4A10", "CD8A",  
                  "CD4", "CCR7", "SELL", "GZMB") 


###ASSIGN CELL TYPES FOR EACH CLUSTER###
#Use files created by filter_and_integrate_data.R:  
  #1. umap_integrated.png  (shows clusters)
  #2. Plots in feature_plots/ (shows expression of marker genes)

new_cluster_ids <- c("CD8+ T naive", "CD4+ T naive", "CD4+ T naive", "CD56 dim NK", "CD14+ Mono", #0
                     "B", "CD4+ Tcm", "CD14+ Mono", "CD8+ Temra","CD14+ Mono",  #5
                     "MAIT", "CD4+ Tcm", "B", "CD8+ Tem", "FCGR3A+ Mono", #10
                     "CD4+ T naive", "CD14+ Mono", "Treg", "CD56 bright NK", "cDC", #15
                     "CD14+ Mono", "Unassigned", "pDC", "CD56 dim NK", "Unassigned", #20
                     "CD14+ Mono", "Unassigned", "Unassigned", "Unassigned", "Unassigned", #25
                     "Unassigned","Unassigned") #30


###READ IN SEURAT OBJECT###
setwd(base_directory)
seu <- readRDS(file=rds_file)
DefaultAssay(seu) <- "integrated"


###LABEL CLUSTERS

#Rename identities using new cluster IDs
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
seu$celltype <- Idents(seu)
my_levels <- unique(new.cluster.ids)
Idents(seu) <- factor(Idents(seu), levels= my_levels)


###MAKE INTEGRATED UMAP WITH NEW CLUSTER LABELS####

#Generate list of colors for new integrated plot  
#Replace "Unassigned" with gray (#D3D3D3)

color_list <- hue_pal()(length(my_levels))
for (i in 1:length(color_list)) {
  if (my_levels[i]=="Unassigned") {
    color_list[i]="#D3D3D3"
  }
}

#Make new labeled UMAP
png("umap_labeled_integrated_original.png", width = 800, height = 600, units = "px")
DimPlot(seu, reduction = "umap", label = TRUE, label.size=4, pt.size = 0.5, cols=color_list, repel=TRUE) + NoLegend() 
dev.off()


###OPTIONAL:  Change the order of the cell types to improve the distribution of colors on UMAP###
  #Otherwise you can comment out these lines (65-72)

#my_levels <- c("CD4+ T naive","CD8+ T naive",  "CD8+ Tem", "CD14+ Mono", "FCGR3A+ Mono", "CD4+ Tcm","CD56 bright NK",
#               "MAIT", "pDC",  "B", "Treg",  "CD56 dim NK", "CD8+ Temra", "Unassigned","cDC")  
#Idents(seu) <- factor(Idents(seu), levels= my_levels)

#png("umap_labeled_integrated_improved.png",width = 800, height = 600, units = "px")
#DimPlot(seu, reduction = "umap", label = TRUE, label.size=4, pt.size = 0.5, cols=color_list, repel=TRUE) + NoLegend() 
#dev.off()


###SAVE LABELED RDS FILE###
saveRDS(seu, file = "integrated_data_processed_labeled.rds")


###DROP UNASSIGNED CELLS###
  #They will still be present in the saved RDS from the previous step, in case you need them in the future
seu <- subset(seu, subset = celltype != "Unassigned")


###OPTIONAL: Rename levels for violin plots###
  #Otherwise you can comment out these lines (85-89)

#my_levels <- c("B", "CD14+ Mono", "FCGR3A+ Mono","cDC", "pDC", 
#               "CD56 bright NK", "CD56 dim NK", "MAIT", "Treg","CD8+ T naive", 
#               "CD8+ Tem", "CD8+ Temra", "CD4+ T naive", "CD4+ Tcm") # "dN T",

#Idents(seu) <- factor(Idents(seu), levels= my_levels)


###MAKE VIOLIN PLOTS SHOWING EXPRESSION OF MARKER GENES IN EACH CELL TYPE###

#Make new directory to hold violin plots
if (!dir.exists("violin_plots")) {
  system(sprintf("mkdir violin_plots"))
}

#Define function
make_violin <- function(gene) {
  png(sprintf("violin_plots/%s.png",gene), width = 1200, height = 600, units = "px")
  print({
    VlnPlot(seu, features = c(gene), pt.size=0.05)
  })
  dev.off()
}

#Run function on every gene
for (gene in marker_genes) {
  make_violin(gene)
}


###MAKE SUMMARY DOT PLOT OF MARKER GENES###
DefaultAssay(seu) <- "RNA"

png(sprintf("dot_plot_integrated.png",gene), width = 1200, height = 600, units = "px")
DotPlot(seu, features = marker_genes, cols = c("blue", "red"), dot.scale = 8) + 
  RotatedAxis()
dev.off()
