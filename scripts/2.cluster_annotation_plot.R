#scRNAseq Making UMAP, clustering, and annotations

#Load libraries ====
library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(future)

#Increase memory for Cell object
options(future.globals.maxSize = 30 * 1024^3)
plan("sequential")

#Load RDS object 
seurat_rds <- readRDS("../data/seurat_filtered.rds")
seurat_rds #24386 features across 154343 samples


#Normalization of data and scaling with SCTransform()
#Normalization of data and scaling with SCTransform() -> Replaces NormalizeData(), FindVariableFeatures(), ScaleData(), calculates Pearson Residuals
seurat_rds <- SCTransform(seurat_rds,
                          vars.to.regress = "mitoRatio",
                          conserve.memory = TRUE,
                          vst.flavor = "v2",
                          verbose = TRUE)

#Dimensionality reduction by PCA and UMAP embedding
#Standard steps in the Seurat workflow for visualization and clustering
#Identifies eigenvectors that explain the most variation in the data
seurat_rds <- RunPCA(seurat_rds, features = VariableFeatures(object = seurat_rds), verbose = FALSE)

#How many PCs to use for clustering with elbow plot
pdf("../plots/elbow_plot.pdf")
ElbowPlot(seurat_rds)
dev.off()

#UMAP
#20-30 PCs starting point
seurat_rds <- FindNeighbors(seurat_rds, dims = 1:30)
seurat_rds <- FindClusters(seurat_rds, resolution = 0.5)
seurat_rds <- RunUMAP(seurat_rds, dims = 1:30)

#Save the UMAP
pdf("../plots/umap_clusters.pdf")
print(DimPlot(seurat_rds, reduction = "umap", label = TRUE))
dev.off()

#Save final processed object
saveRDS(seurat_rds, "../data/seurat_processed.rds")
