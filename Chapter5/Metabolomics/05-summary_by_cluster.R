rm(list = ls())

library(openxlsx)
library(sharp)
library(summarise)

# Parameters
centroid <- TRUE

# Loading the data
covars <- readRDS("Data/Covariates_no_rep.rds")
metab_neg <- readRDS("Data/Metabolites_negative_imputed.rds")
metab_pos <- readRDS("Data/Metabolites_positive_imputed.rds")
covars <- covars[rownames(metab_neg), ]
print(all(rownames(metab_neg) == rownames(metab_pos)))
print(all(rownames(covars) == rownames(metab_pos)))

# Loading the results
stab <- readRDS("Results/2-clustering/consensus_clustering_very_fine_both.rds")

# Restricting to Women
covars <- covars[which(as.character(covars$gender) == "Female"), ]
metab_pos <- metab_pos[rownames(covars), ]
metab_neg <- metab_neg[rownames(covars), ]

# Merging positive and negative mode
metab <- cbind(metab_pos, metab_neg)

# Extracting stable clusters
myclusters <- Clusters(stab)

# Summarising the data by cluster (medoids)
summarised <- GroupingSummary(
  x = metab,
  group = myclusters,
  summary = "medoid"
)
metab_summarised <- summarised$xsummarised
description_summary <- summarised$description
representatives <- colnames(metab_summarised)

# Saving the outputs
write.xlsx(description_summary, "Tables/2-clustering/Description_summary_medoid.xlsx")
saveRDS(metab_summarised, "Data/Metabolites_imputed_summarised_medoid.rds")

# Summarising the data by cluster (centroids)
summarised <- GroupingSummary(
  x = metab,
  group = myclusters,
  summary = "centroid"
)
metab_summarised <- summarised$xsummarised
description_summary <- summarised$description
colnames(metab_summarised) <- gsub("F", "C", colnames(metab_summarised))
proxy <- colnames(metab_summarised)

# Preparing annotation by representative features
annot_by_representative <- representatives
names(annot_by_representative) <- proxy

# Saving the outputs
write.xlsx(description_summary, "Tables/2-clustering/Description_summary_centroid.xlsx")
saveRDS(metab_summarised, "Data/Metabolites_imputed_summarised_centroid.rds")
saveRDS(annot_by_representative, "Data/Annotation_representative_features.rds")
