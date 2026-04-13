#3. DE and UMAP Clustering again

#Load libraries ====
library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(future)

#Load data ====
srt_rds <- readRDS("../data/seurat_processed.rds")

#Set the limit to 100GB (100 * 1024^3 bytes) so R doesn't limit size automatically
options(future.globals.maxSize = 100 * 1024^3)
plan("multisession", workers = 4)
#Run PCA -> Returns Seurat object with PCA calculation stored in reductions slot
#Calculate 40 PCs but really only will use 16 based on Elbow plot
srt_rds <- RunPCA(srt_rds, 
                  assay = "SCT", #Assay PCA is run on is SCTransform
                  npcs = 40, 
                  verbose = TRUE, #Print the top genes associated with high/low loadings for the PCs
                  seed.use= 42, #Set seed to 42, default
                  )

# # Clustering
srt_rds <- FindNeighbors(srt_rds, dims = 1:16) #Set cutoff as 16 PCs
srt_rds <- FindClusters(srt_rds, resolution = 0.4) #34 Clusters at resolution 0.5 was a bit much
srt_rds <- RunUMAP(srt_rds, dims = 1:16)


#Marker Search
all_markers <- FindAllMarkers(srt_rds,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)
#Save
write.csv(all_markers, "../data/cluster_markers_res0.4.csv")
