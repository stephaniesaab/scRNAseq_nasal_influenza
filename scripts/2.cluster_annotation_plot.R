#scRNAseq Making UMAP, clustering, and annotations

#Load libraries ====
library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(future)

#Increase memory for Cell object
options(future.globals.maxSize = 15 * 1024^3)
plan("multisession", workers = 8)

#Load RDS object 
seurat_rds <- readRDS("../data/seurat_filtered.rds")
seurat_rds #24386 features across 154343 samples


#Normalization of data and scaling with SCTransform()
seurat_rds <- SCTransform(seurat_rds, 
                          vars.to.regress = "mitoRatio", 
                          method = "glmGAMPoi",
                          vst.flavor = "v2",
                          verbose = TRUE)

#Dimensionality reduction by PCA and UMAP embedding
#Standard steps in the Seurat workflow for visualization and clustering
#Identifies eigenvectors that explain the most variation in the data
seurat_rds <- RunPCA(seurat_rds, features = VariableFeatures(object = seurat_rds), verbose = FALSE)

#How many PCs to use for clustering with elbow plot
pdf("plots/elbow_plot.pdf")
ElbowPlot(seurat_rds)
dev.off()

#UMAP
#20-30 PCs starting point
seurat_rds <- FindNeighbors(seurat_rds, dims = 1:30)
seurat_rds <- FindClusters(seurat_rds, resolution = 0.5)
seurat_rds <- RunUMAP(seurat_rds, dims = 1:30)
DimPlot(seurat_rds, reduction = "umap", label = TRUE)

# ====================== 
# 
# # Basic log normalization of the data (Could also instead use SCTransform()) and scaling
# pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize")
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# 
# # By default Seurat only scales variable features. Here, we're instead scaling all features (better downstream visualization)
# # The scaling phase is also where we would regress out unwanted sources of variation, e.g. cell cycle stage
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)
# 
# # We run a PCA to produce principal components that can be used to cluster our cells
# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# 
# # Determine the ‘dimensionality’ of the dataset -> Pick some principle components to use for the rest of the script
# ElbowPlot(pbmc)
# 
# # Clustering
# pbmc <- FindNeighbors(pbmc, dims = 1:17) #Set cutoff as 17 PCs
# pbmc <- FindClusters(pbmc, resolution = 0.5)
# 
# # UMAP with clusters displayed
# pbmc <- RunUMAP(pbmc, dims = 1:17)
# DimPlot(pbmc, reduction = "umap", label = TRUE)
# 
# # Try running the above sections again with a different resolution value. How does it change the clustering?
# 
# # Let's try to annotate these clusters
# # find all markers of cluster 1
# cluster1.markers <- FindMarkers(pbmc, ident.1 = 1)
# head(cluster1.markers, n = 5)
# 
# cluster10.markers <- FindMarkers(pbmc, ident.1 = 10)
# head(cluster10.markers, n = 5)
# 
# # In my run, I found that INPP4B, IL7R, CDC14A, ANK3, CAMK4 were markers for cluster 1
# VlnPlot(pbmc, features = c("INPP4B", "IL7R", "CDC14A", "ANK3", "CAMK4"))
# #Genes for 
# #FeaturePlot() is useful to display the features on our UMAP
# FeaturePlot(pbmc, features = c("INPP4B", "IL7R", "CDC14A", "ANK3", "CAMK4"))
# 
# #Let's check some genes known to be associated with CD4+ T Cells
# FeaturePlot(pbmc, features = c("CD4", "CD3D", "CD3E", "CD3G", "TRAC", "CD8A"))
# 
# # Naive or other T cells?
# FeaturePlot(pbmc, features = c("LEF1", "CCR7", "SELL", "IL2RA"))
# 
# 
# # Labeling clusters
# new.cluster.ids <- c("Cluster 0", "CD4+ T cells", "Cluster 2", "Naive CD4+ T cells", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9", "Cluster 10")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# 
# #Manual annotation
# #INPP4B gene is a phosphatase that is a tumor suppressor in cancer cells -> Cluster 0
# #IL7R in lymphoid lineage cells -> T and B cells
# #CDC14A -> auditory hair cells (inner ear)-> cluster 10
# #ANK3 gene is in neurons
# #CAMK4  T lymphocytes, neurons, and male germ cells
# 
# #TRAC -> T cells
# #CD3E -> T cell receptor complex