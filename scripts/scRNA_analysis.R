#scRNAseq_analysis.R

#Load libraries ====
library(dplyr)
library(Seurat)
library(ggplot2)

seurat_obj <- readRDS("seurat_ass4.rds")


#Get Mitochondrial Percentage
# Mouse gene symbols use "mt-", human use "MT-"
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# 2. Calculate Ribosomal Percentage (Optional but helpful for 'Excellent' grade)
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]")

# 3. View the metadata table to see your new columns
head(seurat_obj@meta.data)

metadata <- seurat_obj@meta.data

#Quality Control ====

#Visualize the number of cell counts per sample
#Can get multiple cellular barcodes per hydrogen droplet, depending on protocol
#Can get a higher number of cell barcodes than cells (also account for dying cells)
metadata %>%
  ggplot(aes(x = biosample_id, fill = organ_custom))+
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("NumCells")
#Way more cells in OM and RM organs than LNG
#Probably have some junk cells


#UMI counts per cell
metadata %>%
  ggplot(aes(color = biosample_id, x = nCount_RNA, fill = organ_custom))+
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
#All are 3000 or above -> good

#Genes detected per cell
metadata %>%
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= organ_custom)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = 300)
#Bimodal -> One really high peak for D05OM?

#Complexity
#Ratio of nGenes over nUMI
#If there are many transcripts and low nGenes then likely captured a low number of genes and sequenced transcripts from those repeatedly, could be a specific cell type or artifact or contamination
#Expect generally a novelty score above 0.80 for good quality cells
#Make novelty score

# Add number of genes per UMI for each cell to metadata
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
metadata$log10GenesPerUMI <- seurat_obj$log10GenesPerUMI

metadata %>%
  ggplot(aes(x = log10GenesPerUMI, color = organ_custom, fill = orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
#Almost all are above 0.8 -> Good

#Mitochondrial counts ratio
#Can have mitochondrial contamination from dead or dying cells -> poor quality samples surpass 0.2 mt ratio
#Lowecase "mt" because mice, "MT" for humans
# Compute percent mito ratio
seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^mt-")
metadata$mitoRatio <- seurat_obj$mitoRatio / 100

#Plot mtcounts, removed 89 rows becase out of range for log10 scale
metadata %>%
  ggplot(aes(fill = orig.ident, x = mitoRatio, color = organ_custom)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.2)
#All are below 0.2 -> Good

#Joint filtering ====
#Set thresholds for individual metrics to be as permissive as possible and consider joint effects -> reduce risk of filtering viable cell populations
#Good cells show high #genes and #UMIs (upper right quadrant)
#Bottom right-quadrant: Cells with high UMI but low genes: could be dying cells or low complexity cells
#Darker points have higher mt read fractions in low count cells with few detected genes -> Could be damaged/dying cells

#Two metrics often evaluated together: #UMIs and #Genes
#Plot #genes vs. #UMIs coloured by mt fraction
metadata %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio)) +
  geom_point() +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method = lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~orig.ident)


# Visualize QC metrics as a violin plot
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)


#Filtering ====
#Cell-level filtering thresholds, based on general guidelines and the source paper (Kazer et al.)
#nUMI > 750 < 10,000
#nGene > 500
#log10GenesPerUMI > 0.8
#mitoRatio < 15%
#Filter with Seurat function subset()

#Filter out low quality cells with the thresholds above
filtered_seurat <- subset(x = seurat_obj,
                          subset = (750 < nCount_RNA &
                                      nCount_RNA < 10000) &
                            (nFeature_RNA > 500) &
                            (log10GenesPerUMI > 0.80) &
                            (mitoRatio < 15)) #15%
filtered_seurat  #25129 features across 154343 samples within 1 assay
seurat_obj #25129 features across 156572 samples within 1 assay

#Gene-level filtering ====
#Get rid of genes with zero counts

#Extract counts
counts <- GetAssayData(object = filtered_seurat, layer = "counts")

#Output a logical matrix specifying whether each gene is non-zero count or not
nonzero <- counts > 0

#Filter for genes expressed in 10+ cells
#Sum all TRUE values and return TRUE of > 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

#Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

#Filter seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
filtered_seurat #An object of class Seurat 24386 features across 154343 samples within 1 assay

#Reassess after QC

# Visualize QC metrics as a violin plot
VlnPlot(filtered_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio"), ncol = 3)

#Remove objects to save space
rm(counts, nonzero, keep_genes, filtered_counts)
# Save using RDS
saveRDS(filtered_seurat, file = "../data/seurat_filtered.rds")
cat("Cells before:", ncol(seurat_obj), "\n")
cat("Cells after:", ncol(filtered_seurat), "\n")
cat("Features before:", 25129, "\n")
cat("Cells after:", 24386, "\n")


