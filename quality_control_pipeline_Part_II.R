#Generates quality-control plots separated out by cell type, including violinplots with # features and % mito by celltype, 
#celltype-specific heatmaps with correlation of gene expression between samples, and heatmap showing proportion cell type by sample

library(ggplot2)
library(data.table)
library(Seurat)
library(Matrix)
theme_set(theme_bw())
library(dplyr)
library(cowplot)
library(ggdendro)
library(reshape2)
library(RColorBrewer)
library(gplots)

###Set base directory###
base_dir <- "/lab/page_human_data/emily/single_cell/"

####Read in integrated Seurat object###
setwd(base_dir)
seu <- readRDS(file="integrated_results/integrated_data_processed_labeled.rds")
DefaultAssay(seu) <- "RNA"
seu$celltype <- Idents(seu)

###Drop unassigned cells###
seu <- subset(seu, subset = celltype != "Unassigned")

###Make correlation heatmaps and plots for each cell type###
#https://github.com/satijalab/seurat/issues/1552


#Define correlation heatmap function
make_correlation_heatmap <- function(seu, name) {
  
  Idents(seu) <- seu$orig.ident
  
  #Correlation between samples
  av.exp <- AverageExpression(seu)$RNA
  cor.exp <- as.data.frame(cor(av.exp,method="spearman"))
  cor.exp$x <- rownames(cor.exp)
  print(head(cor.exp))
  cor.df <- tidyr::gather(data = cor.exp, y, correlation, c(levels(unique(Idents(seu)))))
  
  png(sprintf("quality_control/correlation_heatmap_%s.png",name))
  print( {
  ggplot(cor.df, aes(x, y, fill = correlation)) +
    geom_tile() + 
    theme(axis.text.x = element_text(angle = 90, vjust=0,hjust=0)) +
    ylab("") + xlab("")
  })
  dev.off()
}

#Define correlation plot function
make_correlation_plot <- function(seu, name) {
  
  Idents(seu) <- seu$orig.ident
  
  samples <- sample(unique(seu$orig.ident), 2)
  sample_1 <- samples[1]
  sample_2 <- samples[2]
  
  seu_1 <- subset(seu, orig.ident==sample_1)
  seu_2 <- subset(seu, orig.ident==sample_2)
  
  av.exp.1.rank <- rank(AverageExpression(seu_1)$RNA)
  av.exp.2.rank <- rank(AverageExpression(seu_2)$RNA)
  print(av.exp.1.rank[1:5])
  
  df <- data.frame(row.names=rownames(seu_1[["RNA"]]@counts), av.exp.1.rank, av.exp.2.rank)
  colnames(df) <- c("sample_1","sample_2")
  reg<-lm(av.exp.1.rank ~ av.exp.2.rank, data = df)
  r2 <- round(summary(reg)$r.squared,2)
  
  png(sprintf("quality_control/correlation_plot_%s.png", name), width = 480, height = 480, units = "px")
  print( {
    ggplot(data = df, aes(y = sample_2,x=sample_1)) + geom_point() + theme_classic() +
      geom_smooth(method = 'lm',color="red",se=FALSE,linetype=2) + 
      geom_hline(yintercept=0,color="blue",size=1) + geom_vline(xintercept=0,color="blue",size=1) + 
      ggtitle(sprintf("R2 = %s",r2)) +
      ylab(sprintf("Ranks %s",sample_2)) + xlab(sprintf("Ranks %s",sample_1))
  } )
  dev.off()
  
}

###Make correlation heatmap and plot for each celltype###
celltypes <- unique(seu$celltype)

for (i in 1:length(celltypes)) {
  type <- celltypes[i]
  
  if (type!="Unassigned") {
    cells_of_interest <- subset(seu, celltype==type)
    
    make_correlation_heatmap(cells_of_interest, type)
    make_correlation_plot(cells_of_interest, type)
  }
}


###Violin plots by cell type###
Idents(seu) <- seu$celltype

png("quality_control/violin_features_by_celltype.png", width = 1480, height = 480, units = "px")
VlnPlot(seu, features = c("nFeature_RNA"), pt.size=0.01)
dev.off()

png("quality_control/violin_ribosomal_by_celltype.png", width = 1480, height = 480, units = "px")
VlnPlot(seu, features = c("percent.ribo"), pt.size=0.01)
dev.off()

png("quality_control/violin_mito_by_celltype.png", width = 1480, height = 480, units = "px")
VlnPlot(seu, features = c("percent.mt"), pt.size=0.01)
dev.off()



###Clustermap of cell type proportions by individual###

#Get proportions by cell type

table <- table(Idents(seu), seu$orig.ident)
write.table(table, file = "quality_control/celltypes_by_sample.csv", sep = ",")

# Read in data
table <- read.csv(file = "quality_control/celltypes_by_sample.csv", stringsAsFactors = TRUE)

#Convert to proportions for each individual (after dropping unassigned)
table <- t(scale(table, center = FALSE, scale = colSums(table)))

#Get z-scores of proportions across each celltype
#table.scaled <- scale(table, center = TRUE)

png("quality_control/celltype_heatmap.png")
heatmap.2(table,margins = c(8, 8),key=TRUE, col= colorRampPalette(brewer.pal(8, "YlOrRd"))(50), tracecol=NA)
dev.off()

###Plot sex chromosome gene expression###

xist = c("XIST")
x_genes = c("RPS4X","DDX3X","EIF1AX","KDM5C","TXLNG", "USP9X","PRKX", "KDM6A")
y_genes = c("EIF1AY","KDM5D","TXLNGY","RPS4Y1","UTY","DDX3Y","ZFY","NLGN4Y","USP9Y","PRKY","TMSB4Y")

#This adds a metadata column (e.g. "percent platelet"), you should be able to do this for lists of X genes and Y genes
seu[["percent.xist"]] <- PercentageFeatureSet(seu, features = xist, assay="RNA")
seu[["percent.x"]] <- PercentageFeatureSet(seu, features = x_genes, assay="RNA")
seu[["percent.y"]] <- PercentageFeatureSet(seu, features = y_genes, assay="RNA")

png("quality_control/xist.png", width = 1000, height = 800, units="px")
VlnPlot(seu, features = c("percent.xist"), split.by = "Groups", idents = c("B", "CD14+ Mono", "CD4+ T naive"), pt.size = 0)
dev.off()
png("quality_control/xgenes.png", width = 1000, height = 800, units="px")
VlnPlot(seu, features = c("percent.x"), split.by = "Groups", idents = c("B", "CD14+ Mono", "CD4+ T naive"), pt.size = 0)
dev.off()
png("quality_control/ygenes.png", width = 1000, height = 800, units="px")
VlnPlot(seu, features = c("percent.y"), split.by = "Groups", idents = c("B", "CD14+ Mono", "CD4+ T naive"), pt.size = 0)
dev.off()


#UMAPS by individuals
png("quality_control/umap_by_sample.png", width = 1000, height = 1000, units = "px")
DimPlot(seu, reduction = "umap", split.by = "orig.ident", ncol=6 )  + theme(legend.position = 'none')
dev.off()
