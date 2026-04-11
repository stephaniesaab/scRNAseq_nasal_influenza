# ==============================================================================
# Script: get_markers_0.4.R
# Purpose: Generate Full and Top 5 markers for manual annotation (0.4 res)
# ==============================================================================

# 1. LOAD LIBRARIES
library(Seurat)
library(future)
library(dplyr)

# 2. SETUP MULTITHREADING (COMPUTE CANADA OPTIMIZED)
# We use 16 workers to match our Slurm request
plan("sequential")

# Set memory limit to 100GB to prevent "globals" errors with 154k cells
options(future.globals.maxSize = 100 * 1024^3)

message(">>> Sequential active. Starting Marker Discovery...")

# 3. LOAD THE DATA
# Using the object name you provided
srt <- readRDS("../data/seurat_cluster0.4.rds")

# 4. SET IDENTITY
# Ensure the 0.4 clusters are the active identity
# We use the column name identified from your previous metadata check
if ("SCT_snn_res.0.4" %in% colnames(srt@meta.data)) {
    Idents(srt) <- "SCT_snn_res.0.4"
} else {
    # Fallback to seurat_clusters if the specific res column isn't found
    Idents(srt) <- "seurat_clusters"
}

message(paste(">>> Identifying markers for", length(unique(Idents(srt))), "clusters."))

# 5. RUN FINDALLMARKERS
# only.pos = TRUE: We only want upregulated genes (identity markers)
# min.pct = 0.25: Gene must be in at least 25% of cells in the cluster
# logfc.threshold = 0.25: Gene must be at least 1.2x higher than others
all_markers <- FindAllMarkers(
    srt, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    logfc.threshold = 0.25,
    max.cells.per.ident = 1000, verbose = TRUE)

# 6. SAVE THE FULL SET
# This is your "Deep Dive" file
write.csv(all_markers, "../data/cluster_markers_0.4_FULL.csv", row.names = FALSE)
message(">>> Full marker list saved to ../data/cluster_markers_0.4_FULL.csv")

# 7. GENERATE AND SAVE THE TOP 5
# We rank by 'avg_log2FC' to get the most specific markers per cluster
top5_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

write.csv(top5_markers, "../data/cluster_markers_0.4_TOP5.csv", row.names = FALSE)
message(">>> Top 5 marker list saved to ../data/cluster_markers_0.4_TOP5.csv")

# 8. PRINT RESULTS FOR IMMEDIATE VIEWING
# This will show up in your Slurm .out file
print(top5_markers, n = 150)

message(">>> Job Complete. Ready for manual annotation!")
