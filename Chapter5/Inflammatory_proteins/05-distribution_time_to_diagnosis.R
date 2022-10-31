rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter5/")

library(openxlsx)
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

pdf("Working_figures/Inflammatory_proteins/Histogram_time_to_diagnosis.pdf")
par(mar = c(5, 5, 1, 1))
hist(covars$ttd,
  col = "navy", border = "white",
  xlab = "Time to diagnosis (years)",
  ylab = "Count",
  main = "",
  cex.lab = 1.5, las = 1, breaks = 10
)
abline(v = median(covars$ttd, na.rm = TRUE), col = "red", lty = 2, lwd = 2)
legend("topright", lty = 2, lwd = 2, col = "red", legend = "Median", bty = "n", cex = 1.5)
dev.off()
