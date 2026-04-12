#6.DE_GSEA

#Load libraries
library(Seurat)
library(SingleCellExperiment) 
library(SingleR)
library(DESeq2) #For the DE
#library(clusterProfiler) Not in the narval cluster
library(org.Mm.eg.db)        
library(BiocParallel)

#Setup parallelization
plan("multisession", workers = 8) 
options(future.globals.maxSize = 100 * 1024^3) #100GB limit

#Load data
srt <- readRDS("../data/seurat_cluster0.4.rds")

#Rename clusters with cell types annotated manually
new_cluster_ids <- c(
  "0" = "Calb2+ OSN",
  "1" = "Mature OSN",
  "2" = "Fcrls+Trem2+ Mac",
  "3" = "IFN-Stim Basal", #IFN Stim speculated
  "4" = "B-Cell",
  "5" = "Dlg2+ OSN",
  "6" = "Fibroblast",
  "7" = "Neural Progenitor",
  "8" = "IFN-Stim DC", #IFN Stim speculated
  "9" = "Endothelial",
  "10" = "T/NK Cells",
  "11" = "Dlg2+ OSN",
  "12" = "MHC-I Hi Mature",
  "13" = "Ciliated",
  "14" = "Igfbp2+Nrcam+ Basal",
  "15" = "Non-Classical Mono",
  "16" = "Mature OSN",
  "17" = "Gpx6+Ces1a+ Gob/Sec",
  "18" = "Ionocyte",
  "19" = "Serous epithelial",
  "20" = "Immature Neutrophil",
  "21" = "Smooth Muscle Cell",
  "22" = "Sustentacular (OE)",
  "23" = "Scgb-b27+Cck+",
  "24" = "Glandular epithelial",
  "25" = "Nasal Tuft",
  "26" = "Neutrophil Progenitor"
)

#unname() prevents Seurat from looking at the '0', '1', '2' names and trying to match cell barcodes
#Add new cluster IDs to metadata
srt$cell_type <- unname(new_cluster_ids[as.character(srt$SCT_snn_res.0.4)])

#Filter for groups to compare
#GSEA and DE of clusters 3 and 8 compared between sample groups: 
#Only want biosample IDs: Naive_RM, D05_RM, D14_RM for clusters 3 and 8
keep_samples <- c("Naive_RM", "D05_RM", "D14_RM",
                  "Naive_OM", "D05_OM", "D14_OM")

#Filter for cells of interest, this will keep 71101 cells
srt <- subset(srt, subset = biosample_id %in% keep_samples)

#Pseudobulk aggregation
#Sum counts by Sample + CellType + Day + Organ
#Turns thousands of cells into a few pseudobulk samples
pseudo_srt <- AggregateExpression(srt,
                                  assays = "RNA",
                                  return.seurat = TRUE,
                                  group.by = c("cell_type", "organ_custom", 
                                               "time", "biosample_id", "mouse_id"))
ncol(pseudo_srt) #Should have 3x as many columns (3x12)
#Create comparison ID on pseudobulk object
pseudo_srt$cell_time_organ <- paste(pseudo_srt$cell_type,
                                     pseudo_srt$time,
                                     pseudo_srt$organ_custom, sep = "_")
Idents(pseudo_srt) <- "cell_time_organ"

#Running differential expression (DE) loop
#Want to compare clusters 3 and 8 across time points
#Want to compare clusters 5 and 11 (Merged to Mature OSN) in OM
#ct = cell type
#tp = timepoint
#org = organ
#For mouse id naming conventions with paper: RT = RM, ET = OM
comparisons <- list(
  list(ct="IFN-Stim Basal", tp="D05", org="RM"),
  list(ct="IFN-Stim DC",    tp="D05", org="RM"),
  list(ct="Mature OSN", tp="D05", org="OM"),
  list(ct="IFN-Stim Basal", tp="D14", org="RM"),
  list(ct="IFN-Stim DC",    tp="D14", org="RM"),
  list(ct="Mature OSN", tp="D14", org="OM")
)

for (item in comparisons) {
  g1 <- paste0(item$ct, "_", item$tp, "_", item$org)
  g2 <- paste0(item$ct, "_Naive_", item$org)
  
  if (g1 %in% Idents(pseudo_srt) & g2 %in% Idents(pseudo_srt)) {
    
    #For output in cluster when running SLURM
    message(paste(">>> Running DE:", g1, "vs", g2))
    #Run FindMarkers using DESeq2 test
    bulk_de <- Seurat::FindMarkers(object = pseudo_srt, 
                           ident.1 = g1, 
                           ident.2 = g2,
                           logfc.threshold = 0.2, #Default is 0.1, increasing can miss weaker signals but limits tests to genes which show on average x-fold difference between groups of cells
                           test.use = "DESeq2", #negative binomial distribution
                           slot = "counts",#DESeq2 must use raw counts,
                           verbose = TRUE) 
    
    #Save the comparison
    write.csv(bulk_de, paste0("../data/Pseudobulk_DESeq2_", g1, "_vs_Naive.csv"))
   }
}
message(">>> DONE!")
# ====================================
# GSEA
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

#Load DE files from data
de_files <- list.files(path = "../data/", pattern = "^Pseudobulk", full.names = TRUE)

#GSEA Loop
for (file_path in de_files) {
  
  #Get the name for the output file
  base_name <- basename(file_path)
  base_name <- gsub(".csv", "", base_name)
  message(paste(">>> Processing:", base_name))
  
  #Load the DE results using path names
  res_de <- read.csv(file_path, row.names = 1)
  
  #Create the Ranked Gene List
  #Rank by avg_log2FC (Standard for GSEA)
  #Remove any NAs
  res_de <- res_de[!is.na(res_de$avg_log2FC), ]
  
  gene_list <- res_de$avg_log2FC
  names(gene_list) <- rownames(res_de)
  
  #Sort descending
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  #GSEA function
  gsea_res <- gseGO(geneList = gene_list, #Just made
          OrgDb        = org.Mm.eg.db, #Mouse database
          keyType      = 'SYMBOL',
          ont          = "BP", #Biological processes
          pvalueCutoff = 0.05, #Get significant ones
          minGSSize    = 10, #minimal size of each geneSet for analyzing (avoid small pathways that may seem stat sig)
          maxGSSize    = 500, #maximal size of genes annotated for testing, avoid pathways that are too broad (want specifics)
          pAdjustMethod = "BH",
          verbose      = FALSE)
  
  #Save the Results
  if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
    write.csv(as.data.frame(gsea_res), 
              paste0("../data/GSEA_Standalone_", base_name, ".csv"), 
              row.names = FALSE)
    message(paste(">>>Saved GSEA for", base_name))
  } else {
    message(paste("--- No significant pathways found for", base_name))
  }
}
