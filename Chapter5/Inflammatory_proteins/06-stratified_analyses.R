rm(list = ls())

library(openxlsx)
library(sharp)

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

# Checking balanced subsampling by subtype
strata <- ifelse(is.na(covars$subtype), yes = "Control", no = covars$subtype)
ids <- Resample(data = covars$case, family = "binomial", resampling = SubtypeResampling, strata = strata)
names(strata) <- rownames(covars)
prop.table(table(strata))
prop.table(table(strata[ids]))

# Storing full dataset
full_covars <- covars
full_proteins <- proteins
full_strata <- strata

for (stratum in c("never_smokers", "small_ttd", "large_ttd")) {
  print(stratum)

  # Keeping only the stratum
  if (stratum == "never_smokers") {
    covars <- full_covars[which(as.character(full_covars$smoking_status) == "Never"), ]
    covars <- rbind(
      covars[which(covars$case == "0"), ],
      covars[which(tolower(covars$subtype) == "adenocarcinoma"), ]
    )
    proteins <- full_proteins[rownames(covars), ]
    strata <- full_strata[rownames(covars)]
  }
  if (stratum == "small_ttd") {
    ids <- which((full_covars$case == 0) | (full_covars$ttd < median(full_covars$ttd, na.rm = TRUE)))
    covars <- full_covars[ids, ]
    proteins <- full_proteins[rownames(covars), ]
    strata <- full_strata[rownames(covars)]
  }
  if (stratum == "large_ttd") {
    ids <- which((full_covars$case == 0) | (full_covars$ttd >= median(full_covars$ttd, na.rm = TRUE)))
    covars <- full_covars[ids, ]
    proteins <- full_proteins[rownames(covars), ]
    strata <- full_strata[rownames(covars)]
  }

  # Printing numbers
  print(table(covars$case))

  # Associations with lung cancer (multivariate)
  if (stratum != "never_smokers") {
    conf <- c("age.sample", "packyears")
  } else {
    conf <- c("age.sample")
  }
  xdata <- cbind(proteins, covars[, conf, drop = FALSE])
  stab <- VariableSelection(
    xdata = xdata,
    ydata = covars$case,
    penalty.factor = c(rep(1, ncol(proteins)), rep(0, length(conf))),
    family = "binomial",
    resampling = SubtypeResampling,
    strata = strata,
    output_data = TRUE
  )

  # Calibration plot
  par(mar = c(7, 5, 5, 5))
  CalibrationPlot(stab)

  # Explanatory performances
  expl <- Incremental(
    xdata = xdata,
    ydata = covars$case,
    stability = stab,
    n_predictors = ncol(xdata),
    K = 100,
    resampling = SubtypeResampling,
    strata = strata
  )

  # Excluding age
  expl$AUC <- expl$AUC[length(conf):length(expl$AUC)]
  expl$names <- expl$names[length(conf):length(expl$names)]
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
  saveRDS(stab, paste0("Results/Inflammatory_proteins/Multivariate_lung_cancer_", stratum, ".rds"))
  saveRDS(expl, paste0("Results/Inflammatory_proteins/Incremental_performances_lung_cancer_", stratum, ".rds"))
}
