#5.annotating_clusters

#Load libraries
library(SingleR)
library(celldex)
library(BiocParallel)
library(Seurat)
library(SingleCellExperiment)

#Load the mouse reference on login node so it has wifi
ref <- MouseRNAseqData()

#Load seurat object
srt_rds <- readRDS("../data/seurat_processed.rds")

#Make SCE conversion before (required for singleR)
test_sce <- as.SingleCellExperiment(srt_rds)

#Run SingleR on the Seurat object
results <- SingleR(test = test_sce, 
                   ref = ref, 
                   labels = ref$label.main,
		   BPPARAM = MulticoreParam(8)) #'label.main' to get broad cell types and multicore to use 8 cores of CPU

#Add the predicted labels to Seurat object
srt_rds$SingleR_broad <- results$labels

#Find cells with low-confidence (maybe doublets)
#Lower delta means cell looks like two types (two key markers of opposing cell types)
to_prune <- pruneScores(results, min.diff.med = 0.05) 

#Label low-conf cells as 'Ambiguous/Doublet'
srt_rds$SingleR_clean <- srt_rds$SingleR_broad
srt_rds$SingleR_clean[to_prune] <- "Ambiguous_Doublet"

#Write csv to see how SingleR labels compare to my manual cluster numbers
table_compare <- table(srt_rds$seurat_clusters, srt_rds$SingleR_clean)
write.csv(table_compare, "../data/SingleR_vs_Clusters.csv")


#Save the Annotated seurat object
saveRDS(srt_rds, "../data/seurat_auto_anno.rds")

