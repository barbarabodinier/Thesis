rm(list = ls())

library(openxlsx)

source("~/Dropbox/PhD/Thesis/Version2/chapter5/Scripts/functions.R")

# Parameters
all_features <- FALSE
summary_method <- "centroid"

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
covars$packyears <- scale(covars$packyears)
covars$bmi <- scale(covars$bmi)
covars$age.sample <- scale(covars$age.sample)

# Associations with smoking (univariate)
out <- RunLogistic(
  covars = covars,
  conf = c("age.sample", "packyears"),
  outcome = "case",
  predictors = metab
)

# Saving the results
if (all_features) {
  saveRDS(out$summary, "Results/Metabolomics/Univariate_lung_cancer_all_features.rds")
  write.xlsx(out$formatted, paste0("Tables/Metabolomics/Univariate_lung_cancer_all_features.xlsx"),
    rowNames = TRUE, overwrite = TRUE
  )
} else {
  saveRDS(out$summary, "Results/Metabolomics/Univariate_lung_cancer.rds")
  write.xlsx(out$formatted, paste0("Tables/Metabolomics/Univariate_lung_cancer.xlsx"),
    rowNames = TRUE, overwrite = TRUE
  )
}
