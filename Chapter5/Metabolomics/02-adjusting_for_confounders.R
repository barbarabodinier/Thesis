rm(list = ls())

library(sharp)

# Creating folders
filepath <- "Results/2-clustering/"
dir.create("Results", showWarnings = FALSE)
dir.create(filepath, showWarnings = FALSE)

# Loading the data
covars <- readRDS("Data/Covariates_no_rep.rds")
metab_negative <- readRDS("Data/Metabolites_negative_imputed.rds")
metab_positive <- readRDS("Data/Metabolites_positive_imputed.rds")
covars <- covars[rownames(metab_negative), ]
print(all(rownames(metab_negative) == rownames(metab_positive)))
print(all(rownames(covars) == rownames(metab_positive)))

# Kepping only Women
covars <- covars[which(as.character(covars$gender) == "Female"), ]
metab_positive <- metab_positive[rownames(covars), ]
metab_negative <- metab_negative[rownames(covars), ]

# Standardising continuous covariates
covars$packyears <- scale(covars$packyears)
covars$bmi <- scale(covars$bmi)
covars$age.sample <- scale(covars$age.sample)

# Creating data adjusted on smoking and age
for (mode in c("positive", "negative")) {
  print(mode)
  metab <- eval(parse(text = paste0("metab_", mode)))
  metab_adjusted_conf=metab
  pb <- utils::txtProgressBar(style = 3)
  for (j in 1:ncol(metab)){
    metab_adjusted_conf[,j]=residuals(lm(metab[,j]~covars$age.sample+covars$bmi+covars$packyears))
    utils::setTxtProgressBar(pb, j / ncol(metab))
  }
  saveRDS(metab_adjusted_conf, paste0("Data/Metabolites_",mode,"_adjusted_conf.rds"))
}
