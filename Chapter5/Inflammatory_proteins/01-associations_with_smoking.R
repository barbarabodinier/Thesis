rm(list = ls())

library(openxlsx)

source("~/Dropbox/PhD/Thesis/Version2/chapter5/Scripts/functions.R")

# Loading the data
proteins <- readRDS("Data/Original/Proteins_selected_denoised_re.rds")
covars <- readRDS("Data/Original/Covariates_no_rep.rds")

# Excluding participants with missing characteristics
covars <- covars[rownames(proteins), ]
covars <- covars[which(!is.na(covars$bmi) & !is.na(covars$packyears)), ]
proteins <- proteins[rownames(covars), ]
print(nrow(proteins))
print(all(rownames(proteins) == rownames(covars)))

# Keeping only Women
covars <- covars[which(as.character(covars$gender) == "Female"), ]
proteins <- proteins[rownames(covars), ]
print(nrow(proteins))
print(all(rownames(proteins) == rownames(covars)))

# Standardising continuous covariates
covars$packyears <- scale(covars$packyears)
covars$bmi <- scale(covars$bmi)
covars$age.sample <- scale(covars$age.sample)

# Associations with smoking (univariate)
out <- RunLinear(
  covars = covars,
  predictor = "packyears",
  outcome = proteins,
  conf = c("age.sample"),
  multiple_x = FALSE
)

# Saving the results
saveRDS(out$summary, "Results/Inflammatory_proteins/Univariate_packyears.rds")
write.xlsx(out$formatted, paste0("Tables/Inflammatory_proteins/Univariate_packyears.xlsx"),
  rowNames = TRUE, overwrite = TRUE
)
