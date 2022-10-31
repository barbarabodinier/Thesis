rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter5/")

library(openxlsx)
source("~/Dropbox/PhD/Thesis/Version2/chapter5/Scripts/functions.R")

# Loading the data
cpg <- readRDS("Data/Original/imputed_denoised_450K.rds")
covars <- readRDS("Data/Covariates_dna_methylation.rds")

# Excluding participants with missing characteristics
ids <- intersect(rownames(cpg), rownames(covars))
cpg <- cpg[ids, ]
covars <- covars[ids, ]
print(nrow(cpg))
print(all(rownames(cpg) == rownames(covars)))

# # Keeping only Women
# covars <- covars[which(as.character(covars$gender) == "Female"), ]
# cpg <- cpg[rownames(covars), ]
# print(nrow(cpg))
# print(all(rownames(cpg) == rownames(covars)))

# Standardising continuous covariates
covars$packyears <- scale(covars$packyears)
covars$bmi <- scale(covars$bmi)
covars$age.recr <- scale(covars$age.recr)

# Associations with smoking (univariate)
out <- RunLinear(
  covars = covars,
  predictor = "packyears",
  outcome = cpg,
  conf = c("age.recr", "gender"),
  multiple_x = FALSE
)

# Saving the results
saveRDS(out$summary, "Results/DNA_methylation/Univariate_packyears_all.rds")
write.xlsx(out$formatted, paste0("Tables/DNA_methylation/Univariate_packyears_all.xlsx"),
  rowNames = TRUE, overwrite = TRUE
)
