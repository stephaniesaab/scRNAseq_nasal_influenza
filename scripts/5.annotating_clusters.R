#5.annotating_clusters

#Load libraries
library(SingleR)
library(celldex)

#Load the mouse reference
ref <- MouseRNAseqData()

#Load seurat object
srt_rds <- readRDS("../data/seurat_processed.rds")

#Run SingleR on the Seurat object
#Converting to SingleCellExperiment is needed for SingleR
results <- SingleR(test = as.SingleCellExperiment(srt_rds), 
                   ref = ref, 
                   labels = ref$label.main) # 'label.main' to get broad cell types

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

#Plot the automated labels on a UMAP
png("../plots/SingleR_Automated_UMAP.png", width=10, height=8, units="in", res=300)
DimPlot(srt_rds, group.by = "SingleR_clean", label = TRUE, repel = TRUE) + 
  ggtitle("Automated Annotation (SingleR)")
dev.off()