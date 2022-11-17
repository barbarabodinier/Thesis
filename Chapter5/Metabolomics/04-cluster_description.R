rm(list = ls())

library(sharp)

# Creating folders
dir.create("Tables/2-clustering", showWarnings = FALSE)

# Loading the data
metab_neg <- readRDS("Data/Metabolites_negative_imputed.rds")
metab_pos <- readRDS("Data/Metabolites_positive_imputed.rds")

# Loading the results
stab <- readRDS("Results/2-clustering/consensus_clustering_very_fine_both.rds")

# Distribution of cluster sizes
tmp_size <- table(table(Clusters(stab)))
cluster_sizes <- data.frame(
  Size = names(tmp_size),
  Count = as.numeric(tmp_size)
)
write.xlsx(cluster_sizes, "Tables/2-clustering/Distribution_cluster_sizes.xlsx",
  rowNames = FALSE, colNames = TRUE
)

# Membership of bi-modal clusters
myclusters <- sort(Clusters(stab))
clusters_list <- split(names(myclusters), f = myclusters)
bimode_clusters <- sapply(clusters_list, FUN = function(x) {
  npos <- sum(x %in% colnames(metab_pos))
  nneg <- sum(x %in% colnames(metab_neg))
  (npos != 0) & (nneg != 0)
})
mybimodeclusters <- myclusters[myclusters %in% which(bimode_clusters)]
bimode_cluster_membership <- data.frame(
  Cluster_id = as.numeric(mybimodeclusters),
  Feature_name = names(mybimodeclusters),
  Mass = gsub("@.*", "", names(mybimodeclusters)),
  Retention_time = gsub(".*@", "", names(mybimodeclusters)),
  Mode = ifelse(names(mybimodeclusters) %in% colnames(metab_pos),
    yes = "Positive", no = "Negative"
  )
)
write.xlsx(bimode_cluster_membership, "Tables/2-clustering/Bimode_cluster_membership.xlsx",
  rowNames = FALSE, colNames = TRUE
)
