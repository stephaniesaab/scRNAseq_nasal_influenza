#3. DE and UMAP Clustering again

#Load libraries ====
library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(future)

#Load data ====
srt_rds <- readRDS("../data/seurat_processed.rds")

#Run PCA -> Returns Seurat object with PCA calculation stored in reductions slot
srt_rds <- RunPCA(srt_rds, 
                  assay = "SCT", #Assat PCA is run on is SCTransform
                  npcs = 16, #As per elbow plot (Bend around 16PCs)
                  verbose = TRUE, #Print the top genes associated with high/low loadings for the PCs
                  seed.use= 42, #Set seed to 42, default
                  )

# # Clustering
srt_rds <- FindNeighbors(srt_rds, dims = 1:16) #Set cutoff as 17 PCs
srt_rds <- FindClusters(srt_rds, resolution = 0.4) #34 Clusters at resolution 0.5, a bit much

# # UMAP with clusters displayed
srt_rds <- RunUMAP(srt_rds, dims = 1:16)
DimPlot(srt_rds, reduction = "umap", label = TRUE)

#Annotating the clusters
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
# ====================== 

# # The scaling phase is also where we would regress out unwanted sources of variation, e.g. cell cycle stage
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)
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
