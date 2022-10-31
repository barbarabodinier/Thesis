rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter5/")

library(openxlsx)
source("~/Dropbox/PhD/Thesis/Version2/chapter5/Scripts/functions.R")

# Parameters
all_features <- FALSE
summary_method <- "centroid"
smoking_metric <- "ever_smoker"

# Loading the data
if (all_features) {
  covars <- readRDS("Data/Original/Covariates_no_rep.rds")
  metab_neg <- readRDS("/Users/barbara/Dropbox/Metabolomics_lung_cancer/Data/Metabolites_negative_imputed.rds")
  metab_pos <- readRDS("/Users/barbara/Dropbox/Metabolomics_lung_cancer/Data/Metabolites_positive_imputed.rds")
  covars <- covars[rownames(metab_neg), ]
  print(all(rownames(metab_neg) == rownames(metab_pos)))
  print(all(rownames(covars) == rownames(metab_pos)))
  metab <- cbind(metab_pos, metab_neg)
} else {
  metab <- readRDS(paste0("/Users/barbara/Dropbox/Metabolomics_lung_cancer/Data/Metabolites_imputed_summarised_", summary_method, ".rds"))
  covars <- readRDS("Data/Original/Covariates_no_rep.rds")
}

# Excluding participants with missing characteristics
covars <- covars[rownames(metab), ]
covars <- covars[which(!is.na(covars$bmi) & !is.na(covars$packyears)), ]
metab <- metab[rownames(covars), ]
print(nrow(metab))
print(all(rownames(metab) == rownames(covars)))

# Keeping only Women
covars <- covars[which(as.character(covars$gender) == "Female"), ]
metab <- metab[rownames(covars), ]
print(nrow(metab))
print(all(rownames(metab) == rownames(covars)))

# Standardising continuous covariates
if (smoking_metric != "ever_smoker") {
  covars[, smoking_metric] <- scale(covars[, smoking_metric])
}
covars$bmi <- scale(covars$bmi)
covars$age.sample <- scale(covars$age.sample)
metab <- scale(metab)
covars$ever_smoker <- ifelse(as.character(covars$smoking_status) == "Current", yes = 1, no = 0)

# Associations with smoking (univariate)
out <- RunLinear(
  covars = covars,
  predictor = smoking_metric,
  outcome = metab,
  conf = c("age.sample"),
  multiple_x = FALSE
)

# Saving the results
if (all_features) {
  saveRDS(out$summary, paste0("Results/Metabolomics/Univariate_", smoking_metric, "_all_features.rds"))
  write.xlsx(out$formatted, paste0("Tables/Metabolomics/Univariate_", smoking_metric, "_all_features.xlsx"),
    rowNames = TRUE, overwrite = TRUE
  )
} else {
  saveRDS(out$summary, paste0("Results/Metabolomics/Univariate_", smoking_metric, ".rds"))
  write.xlsx(out$formatted, paste0("Tables/Metabolomics/Univariate_", smoking_metric, ".xlsx"),
    rowNames = TRUE, overwrite = TRUE
  )
}
