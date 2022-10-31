rm(list = ls())

library(colorspace)
library(sharp)

# Loading the data
proteins <- readRDS("Data/Filtered_lung_cancer/Proteins.rds")
cpg <- readRDS("Data/Filtered_lung_cancer/DNA_methylation.rds")
metab <- readRDS("Data/Filtered_lung_cancer/Metabolomics.rds")
covars <- readRDS("Data/Filtered_lung_cancer/Covariates.rds")
x <- cbind(proteins, cpg, metab)

# Standardising continuous covariates
covars$packyears <- scale(covars$packyears)
covars$bmi <- scale(covars$bmi)
covars$age.sample <- scale(covars$age.sample)

# Denoising for age, BMI and packyears (in Female only)
x_adjusted_conf <- x
pb <- utils::txtProgressBar(style = 3)
for (j in 1:ncol(x)) {
  x_adjusted_conf[, j] <- residuals(lm(x[, j] ~ covars$age.sample + covars$bmi + covars$packyears))
  utils::setTxtProgressBar(pb, j / ncol(x))
}

# Saving the prepared data
saveRDS(x_adjusted_conf, "Data/Filtered_lung_cancer/OMICs_adjusted_conf.rds")
