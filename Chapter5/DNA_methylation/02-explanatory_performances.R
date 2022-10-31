rm(list = ls())

library(glmnet)
library(sharp)

source("Scripts/functions.R")

# Loading the data
cpg <- readRDS("Data/Original/imputed_denoised_450K.rds")
covars <- readRDS("Data/Covariates_dna_methylation.rds")

# Excluding participants with missing characteristics
ids <- intersect(rownames(cpg), rownames(covars))
cpg <- cpg[ids, ]
covars <- covars[ids, ]
print(nrow(cpg))
print(all(rownames(cpg) == rownames(covars)))

# Standardising continuous covariates
covars$packyears <- scale(covars$packyears)
covars$bmi <- scale(covars$bmi)
covars$age.recr <- scale(covars$age.recr)

# Checking balanced subsampling by subtype
strata <- ifelse(is.na(covars$subtype), yes = "Control", no = as.character(covars$subtype))

# Loading the results
stab <- readRDS(paste0("Results/DNA_methylation/stability_lung_cancer_all.rds"))
CalibrationPlot(stab)

# Explanatory performances
conf <- c("age.recr", "gender", "packyears")
xdata <- cbind(cpg, covars[, conf])
expl <- Incremental(
  xdata = xdata,
  ydata = covars$case_lc,
  stability = stab,
  n_predictors = 100,
  K = 100
)

# Excluding age
expl$AUC <- expl$AUC[3:length(expl$AUC)]
expl$names <- expl$names[3:length(expl$names)]
par(mar = c(8, 5, 1, 1))
IncrementalPlot(expl,
  ylab = "AUC",
  col = c(
    "orange",
    rep("darkred", sum(SelectedVariables(stab))),
    rep("grey", ncol(xdata) - sum(SelectedVariables(stab)))
  ),
  quantiles = c(0.25, 0.75)
)

# Saving the results
saveRDS(expl, "Results/DNA_methylation/Incremental_performances_lung_cancer_all.rds")
