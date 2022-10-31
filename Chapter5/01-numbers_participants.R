rm(list = ls())

# DNA methylation
covars <- readRDS("Data/Original/Covariates_no_rep.rds")
cpg <- readRDS("Data/Original/imputed_denoised_450K.rds")
covars_cpg <- readRDS("Data/Original/covariates_cpg.rds")
covars_cpg <- covars_cpg[which(!is.na(covars_cpg$labnr)), ]
covars_cpg <- covars_cpg[!duplicated(covars_cpg$labnr), ]
cpg <- cpg[which(rownames(cpg) %in% rownames(covars_cpg)), ]
rownames(cpg) <- covars_cpg[rownames(cpg), "labnr"]
table(rownames(cpg) %in% rownames(covars))

ids_cpg <- rownames(cpg)[rownames(cpg) %in% rownames(covars)]
table(covars[ids_cpg, "cohort"], covars[ids_cpg, "case"])

table(covars[ids_cpg, "cohort"])
table(is.na(covars$packyears) | is.na(covars$bmi))
table(is.na(covars$packyears) | is.na(covars$bmi), rownames(covars) %in% ids_cpg)
table(is.na(covars$smok_duration) | is.na(covars$smok_intensity) | is.na(covars$bmi), rownames(covars) %in% ids_cpg)
table(is.na(covars$smok_duration) | is.na(covars$smok_intensity) | is.na(covars$bmi), rownames(covars) %in% ids_cpg, covars$cohort, covars$case)

table(is.na(covars$packyears) | is.na(covars$bmi), covars$cohort)

rownames(covars_cpg) <- covars_cpg$labnr
covars_cpg <- covars_cpg[which(rownames(covars_cpg) %in% rownames(covars)), ]

# Proteins
proteins <- readRDS("Data/Original/Proteins_selected_denoised_re.rds")
covars <- readRDS("Data/Original/Covariates_no_rep.rds")

covars <- covars[rownames(proteins), ]
table(is.na(covars$packyears) | is.na(covars$bmi), covars$cohort, covars$case)

# Metabolomics
metab <- readRDS("Data/Original/Metabolites_negative_imputed.rds")
covars <- readRDS("Data/Original/Covariates_no_rep.rds")

covars <- covars[rownames(metab), ]
table(is.na(covars$packyears) | is.na(covars$bmi), covars$cohort, covars$case)

# OMICs integration
covars <- readRDS("Data/Original/Covariates_no_rep.rds")
ids <- intersect(intersect(rownames(covars_cpg), rownames(proteins)), rownames(metab))
covars <- covars[ids, ]
table(is.na(covars$packyears) | is.na(covars$bmi), covars$cohort, covars$case)
