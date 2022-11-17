rm(list = ls())

library(pheatmap)
library(openxlsx)
library(sharp)
library(summarise)

# Parameters
summary_method <- "centroid"

# Loading the raw data
covars <- readRDS("Data/Covariates_no_rep.rds")
metab_neg <- readRDS("Data/Metabolites_negative_imputed.rds")
metab_pos <- readRDS("Data/Metabolites_positive_imputed.rds")
covars <- covars[rownames(metab_neg), ]
annot <- readRDS("Data/Annotation_representative_features.rds")
print(all(rownames(metab_neg) == rownames(metab_pos)))
print(all(rownames(covars) == rownames(metab_pos)))

# Restricting to Women
covars <- covars[which(as.character(covars$gender) == "Female"), ]
metab_pos <- metab_pos[rownames(covars), ]
metab_neg <- metab_neg[rownames(covars), ]

# Merging positive and negative mode
metab <- cbind(metab_pos, metab_neg)

# Loading the results
stab <- readRDS("Results/2-clustering/consensus_clustering_very_fine_both.rds")

# Calibrated consensus matrix
coprop <- ConsensusMatrix(stability = stab)

# Extracting stable clusters from hierarchical clustering
shclust <- stats::hclust(stats::as.dist(1 - coprop), method = "complete")
ordered_labels <- shclust$labels[shclust$order]

# Heatmap of correlations ordered by consensus clustering
metab <- metab[, ordered_labels]
pheatmap(cor(metab),
  cluster_rows = FALSE, cluster_cols = FALSE,
  # cluster_rows = shclust, cluster_cols = shclust,
  show_rownames = FALSE, show_colnames = FALSE,
  border = NA, breaks = seq(-1, 1, length = 100),
  col = colorRampPalette(c("navy", "white", "darkred"))(100),
  filename = "Figures/2-clustering/Heatmap_correlations_before.png"
)

# Loading the summarised data
metab_summarised <- readRDS(paste0("Data/Metabolites_imputed_summarised_", summary_method, ".rds"))

# Heatmap of correlations ordered by consensus clustering
colnames(metab_summarised) <- annot[colnames(metab_summarised)]
ordered_labels <- ordered_labels[ordered_labels %in% colnames(metab_summarised)]
metab_summarised <- metab_summarised[, ordered_labels]
pheatmap(cor(metab_summarised),
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  border = NA, breaks = seq(-1, 1, length = 100),
  col = colorRampPalette(c("navy", "white", "darkred"))(100),
  filename = paste0("Figures/2-clustering/Heatmap_correlations_after_", summary_method, ".png")
)
