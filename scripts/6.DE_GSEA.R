#6.DE_GSEA

#Load libraries
library(Seurat)
library(SingleCellExperiment) 
library(SingleR)
library(DESeq2) #For the DE
#library(clusterProfiler) Not in the narval cluster
library(org.Mm.eg.db)        
library(BiocParallel)

#DE ==================
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
#saveRDS(srt, file = "../data/seurat_cluster0.4_annotated.rds")

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
    message(paste("Running DE:", g1, "vs", g2))
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
message("done")

# GSEA ====================
#Done on own PC and not SLURM 
library(clusterProfiler) #Not available on the r-bundle-bioconductor on Narval Cluster
library(org.Mm.eg.db)
library(dplyr)

#Load DE files from data
de_files <- list.files(path = "../data/", pattern = "^Pseudobulk", full.names = TRUE)

#GSEA Loop
for (file_path in de_files) {
  
  #Get the name for the output file
  base_name <- basename(file_path)
  base_name <- gsub(".csv", "", base_name)
  message(paste("Processing:", base_name))
  
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
    message(paste("Saved GSEA for", base_name))
  } else {
    message(paste("No sig", base_name))
  }
}

#Feature plots ====================
library(patchwork)
library(Seurat)
library(dplyr)
library(ggplot2)

#Load annotated Seurat object (with manually annotated cell clusters)
srt <- readRDS("../data/seurat_cluster0.4_annotated.rds")

#Define gene lists of interest for each cluster of interest -> Based on top 5 markers
genes_c8 <- c("Cd209a", "Ifi205", "Phf11a", "Tnfsf9", "Ms4a4c") #Cluster IFN-Stim DC
genes_c3 <- c("Acaa1b", "Serpinb5", "Krt15", "Defb1", "Anxa8") #Cluster IFN-Stim Basal
genes_osn <- c("Dlg2", "S100a5", "Rgs7", "Pcp4l1", "Kcnmb3") #Cluster Mature OSN

#Create FeaturePlots
#'stack = TRUE' -> separate plots
plot_c8 <- FeaturePlot(srt, features = genes_c8, ncol = 3, pt.size = 0.5) + 
  plot_annotation(title = "IFN Stimulated DC Markers")

plot_c3 <- FeaturePlot(srt, features = genes_c3, ncol = 3, pt.size = 0.5) + 
  plot_annotation(title = "IFN Stimulated Basal Markers")

plot_osn <- FeaturePlot(srt, features = genes_osn, ncol = 3, pt.size = 0.5) + 
  plot_annotation(title = "Mature OSN Markers")

#Save the plots
ggsave("../plots/FeaturePlot_Cluster8.png", plot_c8, width = 12, height = 8)
ggsave("../plots/FeaturePlot_Cluster3.png", plot_c3, width = 12, height = 8)
ggsave("../plots/FeaturePlot_OSNs.png", plot_osn, width = 12, height = 8)

#GSEA Visuals =================
#Function to make a clean GSEA DotPlot from the 6 CSV
#Ran on interactive Salloc session (not slurm job) because quick
plot_my_gsea <- function(csv_path, title_name) {
  df <- read.csv(csv_path)
  
  #Filter for the top 10 upregulated and top 10 downregulated pathways
  top_up <- df %>% 
    filter(NES > 0) %>% #NES > 0 = upregulated
    arrange(p.adjust, desc(NES)) %>% 
    head(5) #Only want top 5 genes, even if there are ties, for visual clarity
  
  top_down <- df %>% 
    filter(NES < 0) %>% #NES < 0 = downregulated
    arrange(p.adjust, NES) %>% #most negative first
    head(5)
    
  plot_df <- rbind(top_up, top_down)
  
  #Plot
  p <- ggplot(plot_df, aes(x = NES, 
                      y = reorder(Description, NES), 
                      color = NES, #Colour by NES
                      size = -log10(p.adjust))) + #P-values are all very small
    geom_point() +
    scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) + #Set colour gradient for NES
    theme_minimal() +
    labs(title = title_name, #Given in function
         y = "Gene Ontology (BP)", #From GSEA function
         x = "Normalized Enrichment Score (NES)",
         color = "Enrichment (NES)",
         size = "-log10(p.adj)") +
    geom_vline(xintercept = 0, linetype = "dashed") #Separate up and down reg
  return(p)
}
#Colour by p-value
plot_my_gsea2 <- function(csv_path, title_name) {
  df <- read.csv(csv_path)
  
  #Filter for the top 10 upregulated and top 10 downregulated pathways
  top_up <- df %>% 
    filter(NES > 0) %>% #NES > 0 = upregulated
    arrange(p.adjust, desc(NES)) %>% 
    head(5) #Only want top 5 genes, even if there are ties, for visual clarity
  
  top_down <- df %>% 
    filter(NES < 0) %>% #NES < 0 = downregulated
    arrange(p.adjust, NES) %>% #most negative first
    head(5)
  
  plot_df <- rbind(top_up, top_down)
  
  #Plot
  p <- ggplot(plot_df, aes(x = NES, 
                           y = reorder(Description, NES), 
                           color = p.adjust, #Colour by p-value
                           size = setSize)) + #Size by number of genes in data found in GO pathway
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") + #Set colour gradient for NES
    theme_minimal() +
    labs(title = title_name, #Given in function
         y = "Gene Ontology (BP)", #From GSEA function
         x = "Normalized Enrichment Score (NES)",
         color = "p.adjust",
         size = "Number of genes") +
    geom_vline(xintercept = 0, linetype = "dashed") #Separate up and down reg
  return(p)
}
#Key comparisons
p1 <- plot_my_gsea("../data/GSEA_Standalone_Pseudobulk_DESeq2_IFN-Stim Basal_D05_RM_vs_Naive.csv", "GSEA: IFN-Stim Basal (D05 vs Naive)")
p2 <- plot_my_gsea("../data/GSEA_Standalone_Pseudobulk_DESeq2_IFN-Stim Basal_D14_RM_vs_Naive.csv", "GSEA: IFN-Stim Basal (D14 vs Naive)")
p3 <- plot_my_gsea("../data/GSEA_Standalone_Pseudobulk_DESeq2_IFN-Stim DC_D05_RM_vs_Naive.csv", "GSEA: IFN-Stim DC (D05 vs Naive)")
p4 <- plot_my_gsea2("../data/GSEA_Standalone_Pseudobulk_DESeq2_IFN-Stim DC_D14_RM_vs_Naive.csv", "GSEA: IFN-Stim DC (D14 vs Naive)")

ggsave("../plots/GSEA_DotPlot_Basal3_D05.png", p1, width = 10, height = 7)
ggsave("../plots/GSEA_DotPlot_Basal3_D14.png", p2, width = 10, height = 7)
ggsave("../plots/GSEA_DotPlot_DC8_D05.png", p3, width = 10, height = 7)
ggsave("../plots/GSEA_DotPlot_DC8_D14.png", p4, width = 10, height = 7)
