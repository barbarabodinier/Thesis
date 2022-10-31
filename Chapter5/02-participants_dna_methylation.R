rm(list = ls())

# DNA methylation
covars <- readRDS("Data/Original/Covariates_no_rep.rds")
covars_cpg <- readRDS("Data/Original/covariates_cpg.rds")
table(rownames(covars_cpg) == covars_cpg$ID)

# Excluding participants with missing packyears or BMI
covars <- covars[which((!is.na(covars$packyears)) & (!is.na(covars$bmi))), ]

# Identifying controls from other cohorts
ids_bc <- which(covars_cpg$lc_bc == "bc")
covars_bc <- covars_cpg[ids_bc, ]
dim(covars_bc)

# Excluding participants with missing packyears or BMI
covars_bc <- covars_bc[which((!is.na(covars_bc$packyears)) & (!is.na(covars_bc$bmi))), ]
dim(covars_bc)
table(covars_bc$cohort, useNA = "always")

# Identifying participants from nested case-control studies of lung cancer
covars_cpg <- covars_cpg[which(!is.na(covars_cpg$labnr)), ]
covars_cpg <- covars_cpg[!duplicated(covars_cpg$labnr), ]
ids_lc <- which(covars_cpg$labnr %in% rownames(covars))
covars_lc <- covars_cpg[ids_lc, ]
table(covars_lc$cohort, useNA = "always")

# Concatenating data from the two studies
final_covars <- rbind(covars_lc, covars_bc)
table(final_covars$cohort, useNA = "always")
table(final_covars$case_lc, useNA = "always")
table(final_covars$cohort, final_covars$case_lc, useNA = "always")
