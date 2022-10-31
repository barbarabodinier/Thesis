rm(list = ls())

library(openxlsx)
library(glmnet)
library(sharp)

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

# Storing full dataset
full_proteins <- proteins
full_covars <- covars

for (subtype in tolower(unique(covars$subtype[!is.na(covars$subtype)]))) {
  print(subtype)

  # Keeping controls and subtype-specific cases
  covars <- full_covars
  proteins <- full_proteins
  covars <- rbind(
    covars[which(covars$case == "0"), ],
    covars[which(tolower(covars$subtype) == subtype), ]
  )
  proteins <- proteins[rownames(covars), ]

  # Associations with lung cancer (multivariate)
  conf <- c("age.sample", "packyears")
  xdata <- cbind(proteins, covars[, conf])
  stab <- VariableSelection(
    xdata = xdata,
    ydata = covars$case,
    penalty.factor = c(rep(1, ncol(proteins)), rep(0, length(conf))),
    family = "binomial",
    output_data = TRUE
  )
  print(which(SelectedVariables(stab) == 1))

  # Calibration plot
  par(mar = c(7, 5, 5, 5))
  CalibrationPlot(stab)

  # Explanatory performances
  expl <- Incremental(
    xdata = xdata,
    ydata = covars$case,
    stability = stab,
    n_predictors = ncol(xdata),
    K = 100
  )

  # Excluding age and BMI
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
  saveRDS(stab, paste0(
    "Results/Inflammatory_proteins/Multivariate_",
    gsub(" ", "_", subtype),
    ".rds"
  ))
  saveRDS(expl, paste0(
    "Results/Inflammatory_proteins/Incremental_performances_",
    gsub(" ", "_", subtype),
    ".rds"
  ))
}
