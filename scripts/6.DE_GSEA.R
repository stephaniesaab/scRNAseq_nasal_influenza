#6.DE_GSEA

#Load libraries
library(Seurat)
library(SingleCellExperiment) 
library(SingleR)
library(DESeq2) #For the DE
library(clusterProfiler)   
library(org.Mm.eg.db)        
library(BiocParallel)

#Setup parallelization
plan("multisession", workers = 16) 
options(future.globals.maxSize = 100 * 1024^3) #100GB limit

#Load data
srt <- readRDS("../data/seurat_auto_anno.rds")

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

# Apply and save to metadata
srt$cell_type <- cell_type_map[as.character(srt$seurat_clusters)]

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
                                               "time", "biosample_id"))

#Create comparison ID on pseudobulk object
pseudo_srt$cell_time_organ <- paste0(pseudo_srt$cell_type, "_", 
                                     pseudo_srt$dpi, "_", 
                                     pseudo_srt$organ)
Idents(pseudo_srt) <- "cell_time_organ"

#Running differential expression (DE) loop
#Want to compare clusters 3 and 8 across time points
#Want to compare clusters 5 and 11 (Merged to Mature OSN) in OM
#ct = cell type
#tp = timepoint
#org = organ
comparisons <- list(
  list(ct="Basal_3", tp="D05", org="RM"),
  list(ct="DC_8",    tp="D05", org="RM"),
  list(ct="OSN_Merged", tp="D05", org="OM"),
  list(ct="Basal_3", tp="D14", org="RM"),
  list(ct="DC_8",    tp="D14", org="RM"),
  list(ct="OSN_Merged", tp="D14", org="OM")
)

for (item in comparisons) {
  g1 <- paste0(item$ct, "_", item$tp, "_", item$org)
  g2 <- paste0(item$ct, "_Naive_", item$org)
  
  if (g1 %in% Idents(pseudo_srt) & g2 %in% Idents(pseudo_srt)) {
    
    #Run FindMarkers using DESeq2 test
    bulk_de <- FindMarkers(object = pseudo_srt, 
                           ident.1 = g1, 
                           ident.2 = g2,
                           test.use = "DESeq2",
                           slot = "counts") #DESeq2 must use raw counts
    
    #Save the comparison
    write.csv(bulk_de, paste0("../data/Pseudobulk_DESeq2_", g1, "_vs_Naive.csv"))
    
    #GSEA (Ranked by log2FoldChange from DESeq2)
    gene_list <- sort(setNames(bulk_de$avg_log2FC, 
                               rownames(bulk_de)), 
                      decreasing = TRUE)
    gsea_res <- gseGO(geneList = gene_list, 
                      OrgDb = org.Mm.eg.db, 
                      keyType = 'SYMBOL', 
                      ont = "BP")
    
    if(nrow(as.data.frame(gsea_res)) > 0) {
      write.csv(as.data.frame(gsea_res), paste0("../data/GSEA_Pseudobulk_", g1, ".csv"))
    }
  }
}


# 
# # We set logfc.threshold to 0 to get ALL genes for GSEA ranking
# Idents(srt) <- "cell_condition"
# de_results <- FindMarkers(srt, 
#                           ident.1 = "Basal_Influenza", 
#                           ident.2 = "Basal_Mock",
#                           logfc.threshold = 0, # Required for GSEA
#                           min.pct = 0.1)
# 
# write.csv(de_results, "../data/DE_Basal_Full_List.csv")
# 
# # 4. PREPARE GENE LIST FOR GSEA
# # We rank genes by their avg_log2FC
# gene_list <- de_results$avg_log2FC
# names(gene_list) <- rownames(de_results)
# gene_list <- sort(gene_list, decreasing = TRUE)
# 
# # 5. RUN GSEA
# # This tells you which biological pathways are "enriched" at the top or bottom of your list
# gsea_res <- gseGO(geneList      = gene_list,
#                   OrgDb         = org.Mm.eg.db,
#                   keyType       = 'SYMBOL',
#                   ont           = "BP", # Biological Process
#                   minGSSize     = 10,
#                   maxGSSize     = 500,
#                   pvalueCutoff  = 0.05,
#                   verbose       = FALSE)
# 
# # 6. SAVE RESULTS & PLOTS
# write.csv(as.data.frame(gsea_res), "../data/GSEA_Basal_Pathways.csv")
# 
# # Generate the "Running Sum" plot for the top pathway
# png("../plots/GSEA_Top_Pathway.png", width=10, height=8, units="in", res=300)
# gseaplot(gsea_res, geneSetID = 1, by = "runningScore", title = gsea_res$Description[1])
# dev.off()
# 
# # Save the final annotated UMAP
# png("../plots/Final_Annotated_UMAP.png", width=10, height=8, units="in", res=300)
# DimPlot(srt_rds, reduction = "umap", label = TRUE, repel = TRUE) + 
#   ggtitle("Annotated Nasal Mucosa Clusters")
# dev.off()
# 
# # Create a combined identity of CellType + Condition
# srt_rds$cell_condition <- paste0(Idents(srt_rds), "_", srt_rds$condition)
# Idents(srt_rds) <- "cell_condition"
# 
# # Run DE: Influenza vs. Mock specifically in Basal cells (Cluster 3)
# de_results <- FindMarkers(srt_rds, 
#                           ident.1 = "IFN-Stim Basal_Influenza", 
#                           ident.2 = "IFN-Stim Basal_Mock")
# 
# write.csv(de_results, "../data/DE_Basal_Infection_vs_Mock.csv")
# 
# library(clusterProfiler)
# library(org.Mm.eg.db)
# 
# #USE SCTRANSFORM ASSAY FOR THE UMAP
# 
# # Get the genes that are significantly UP in the infected group
# up_genes <- rownames(de_results[de_results$p_val_adj < 0.05 & de_results$avg_log2FC > 1, ])
# 
# # Run the enrichment
# ego <- enrichGO(gene = up_genes, 
#                 OrgDb = org.Mm.eg.db, 
#                 keyType = 'SYMBOL', 
#                 ont = "BP")
# 
# # Plot it
# png("../plots/Pathway_Enrichment_Basal.png", width=10, height=8, units="in", res=300)
# dotplot(ego, showCategory=15) + ggtitle("Pathways Up in Infected Basal Cells")
# dev.off()
# 
# # Calculate proportions
# comp_table <- table(Idents(srt_rds), srt_rds$condition)
# prop_table <- prop.table(comp_table, margin = 2)
# 
# # Save the barplot
# png("../plots/Cell_Composition_Change.png", width=8, height=6, units="in", res=300)
# barplot(prop_table, col=rainbow(nrow(prop_table)), 
#         legend = rownames(prop_table), 
#         args.legend = list(x = "topright", bty = "n", inset=c(-0.2, 0)),
#         main="Tissue Shift: Mock vs Influenza")
# dev.off()
# 
# # 1. SETUP PARALLELIZATION
# # This allows Seurat to use the multiple CPUs you request in Slurm
# plan("multisession", workers = 16) 
# options(future.globals.maxSize = 100 * 1024^3) # 100GB limit
# 
# # 2. LOAD DATA
# srt <- readRDS("../data/seurat_final_annotated.rds")
# 
# # 3. RUN DIFFERENTIAL EXPRESSION (DE)
# # Example: Influenza vs Mock in Basal Cells
# # We set logfc.threshold to 0 to get ALL genes for GSEA ranking
# Idents(srt) <- "cell_condition"
# de_results <- FindMarkers(srt, 
#                           ident.1 = "Basal_Influenza", 
#                           ident.2 = "Basal_Mock",
#                           logfc.threshold = 0, # Required for GSEA
#                           min.pct = 0.1)
# 
# write.csv(de_results, "../data/DE_Basal_Full_List.csv")
# 
# # 4. PREPARE GENE LIST FOR GSEA
# # We rank genes by their avg_log2FC
# gene_list <- de_results$avg_log2FC
# names(gene_list) <- rownames(de_results)
# gene_list <- sort(gene_list, decreasing = TRUE)
# 
# # 5. RUN GSEA
# # This tells you which biological pathways are "enriched" at the top or bottom of your list
# gsea_res <- gseGO(geneList      = gene_list,
#                   OrgDb         = org.Mm.eg.db,
#                   keyType       = 'SYMBOL',
#                   ont           = "BP", # Biological Process
#                   minGSSize     = 10,
#                   maxGSSize     = 500,
#                   pvalueCutoff  = 0.05,
#                   verbose       = FALSE)
# 
# # 6. SAVE RESULTS & PLOTS
# write.csv(as.data.frame(gsea_res), "../data/GSEA_Basal_Pathways.csv")
# 
# # Generate the "Running Sum" plot for the top pathway
# png("../plots/GSEA_Top_Pathway.png", width=10, height=8, units="in", res=300)
# gseaplot(gsea_res, geneSetID = 1, by = "runningScore", title = gsea_res$Description[1])