rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter5/")

# Loading the data
proteins <- readRDS("Data/Original/Proteins_selected_denoised_re.rds")
metab <- readRDS("Data/Original/Metabolites_positive_imputed.rds")
cpg <- readRDS("Data/Original/imputed_denoised_450K.rds")
covars <- readRDS("Data/Original/Covariates_no_rep.rds")
covars_cpg <- readRDS("Data/Original/covariates_cpg.rds")

# Changing participant IDs in DNA methylation
covars_cpg <- covars_cpg[which(!is.na(covars_cpg$labnr)), ]
covars_cpg <- covars_cpg[!duplicated(covars_cpg$labnr), ]
cpg <- cpg[which(rownames(cpg) %in% rownames(covars_cpg)), ]
rownames(cpg) <- covars_cpg[rownames(cpg), "labnr"]

# Identifying participants with all measurements
ids <- intersect(rownames(proteins), rownames(metab))
ids <- intersect(ids, rownames(cpg))
covars <- covars[ids, ]
print(table(is.na(covars$bmi)))
table(is.na(covars$packyears))
table(covars[, c("case", "cohort")])
