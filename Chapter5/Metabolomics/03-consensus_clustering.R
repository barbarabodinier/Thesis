rm(list = ls())

library(sharp)

# Creating folders
filepath <- "Results/2-clustering/"
dir.create("Results", showWarnings = FALSE)
dir.create(filepath, showWarnings = FALSE)

# Parameters
nc_min <- 500
nc_max <- 1200
broad_step <- 50
intermediate_step <- 10
fine_step <- 2
very_fine_step <- 1
mode <- "both"

# Loading the data
covars <- readRDS("Data/Covariates_no_rep.rds")
metab_neg <- readRDS("Data/Metabolites_negative_adjusted_conf.rds")
metab_pos <- readRDS("Data/Metabolites_positive_adjusted_conf.rds")
covars <- covars[rownames(metab_neg), ]
print(all(rownames(metab_neg) == rownames(metab_pos)))
print(all(rownames(covars) == rownames(metab_pos)))

# Restricting to Women
covars <- covars[which(as.character(covars$gender) == "Female"), ]
metab_pos <- metab_pos[rownames(covars), ]
metab_neg <- metab_neg[rownames(covars), ]

# Mode-specific clustering
metab <- cbind(metab_pos, metab_neg)

# Broad calibration
print(system.time({
  stab_broad <- Clustering(
    xdata = metab,
    nc = seq(nc_min, nc_max, by = broad_step),
    row = FALSE
  )
}))
CalibrationPlot(stab_broad)

# Saving the output
saveRDS(stab_broad, paste0(filepath, "consensus_clustering_broad_", mode, ".rds"))
argmax <- Argmax(stab_broad)[1]
rm(stab_broad)

# Intermediate calibration
print(system.time({
  stab_intermediate <- Clustering(
    xdata = metab,
    nc = seq(argmax - broad_step,
      argmax + broad_step,
      by = intermediate_step
    ),
    row = FALSE
  )
}))
CalibrationPlot(stab_intermediate)

# Saving the output
saveRDS(stab_intermediate, paste0(filepath, "consensus_clustering_intermediate_", mode, ".rds"))
argmax <- Argmax(stab_intermediate)[1]
rm(stab_intermediate)

# Fine calibration
print(system.time({
  stab_fine <- Clustering(
    xdata = metab,
    nc = seq(argmax - intermediate_step,
      argmax + intermediate_step,
      by = fine_step
    ),
    row = FALSE
  )
}))
CalibrationPlot(stab_fine)

# Saving the output
saveRDS(stab_fine, paste0(filepath, "consensus_clustering_fine_", mode, ".rds"))
argmax <- Argmax(stab_fine)[1]
rm(stab_fine)

# Very fine calibration
print(system.time({
  stab_very_fine <- Clustering(
    xdata = metab,
    nc = seq(argmax - fine_step,
      argmax + fine_step,
      by = very_fine_step
    ),
    row = FALSE
  )
}))
CalibrationPlot(stab_very_fine)

# Saving the output
saveRDS(stab_very_fine, paste0(filepath, "consensus_clustering_very_fine_", mode, ".rds"))
rm(stab_very_fine)
