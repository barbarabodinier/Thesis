rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter5/")

library(glmnet)
library(sharp)
source("~/Dropbox/PhD/Thesis/Version2/chapter5/Scripts/functions.R")

# Parameters
summary_method <- "centroid"

# Loading the data
metab <- readRDS(paste0("/Users/barbara/Dropbox/Metabolomics_lung_cancer/Data/Metabolites_imputed_summarised_", summary_method, ".rds"))
covars <- readRDS("Data/Original/Covariates_no_rep.rds")

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

# Checking balanced subsampling by subtype
strata <- ifelse(is.na(covars$subtype), yes = "Control", no = covars$subtype)
ids <- Resample(data = covars$case, family = "binomial", resampling = SubtypeResampling, strata = strata)
prop.table(table(strata))
prop.table(table(strata[ids]))

# Associations with lung cancer (multivariate)
conf <- c("age.sample", "packyears")
xdata <- cbind(metab, covars[, conf])
stab <- VariableSelection(
  xdata = xdata,
  ydata = covars$case,
  Lambda = LambdaSequence(lmax = 0.08, lmin = 1e-10, cardinal = 100),
  penalty.factor = c(rep(1, ncol(metab)), rep(0, length(conf))),
  family = "binomial",
  K = 1000,
  resampling = SubtypeResampling,
  strata = strata
)

# Calibration plot
par(mar = c(7, 5, 5, 5))
CalibrationPlot(stab)
table(SelectedVariables(stab))

# Explanatory performances
expl <- Incremental(
  xdata = xdata,
  ydata = covars$case,
  stability = stab,
  n_predictors = 150,
  K = 100,
  resampling = SubtypeResampling,
  strata = strata
)

# Excluding age
expl$AUC <- expl$AUC[2:length(expl$AUC)]
expl$names <- expl$names[2:length(expl$names)]
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
saveRDS(stab, "Results/Metabolomics/Multivariate_lung_cancer.rds")
saveRDS(expl, "Results/Metabolomics/Incremental_performances_lung_cancer.rds")
