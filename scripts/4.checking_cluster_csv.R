#4.Checking_cluster_csv
#Checking the large CSV of all the marker genes for each cluster, arrange by top
#markers for each cluster, make a smaller CSV to look at and pick clusters and 
#important genes manually

#Load libraries
library(dplyr)

#Load data
markers <- read.csv("../data/cluster_markers_res0.4.csv")

#Filter and arrange by descending avg_log2FC
top10_markers <- markers %>% 
  group_by(cluster) %>% #Want to get top markers for each cluster
  filter(p_val_adj < 0.05) %>%  #Get only significant markers
  slice_max(n = 10, order_by, avg_log2FC) #Get top 10 markers for each cluster

#Print it in the console
print(top10_markers)

#Save it to view a smaller sliced CSV (outside of compute canada cluster)
write.csv(top10_markers, "../data/top10_markers_clusters.csv")

#To open an interactive session in Narval:
#salloc --time=01:00:00 --ntasks=1 --cpus-per-task=1 --mem=4G --account=def-itobias
