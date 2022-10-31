rm(list = ls())

library(sharp)

# Loading the data
covars_cpg <- readRDS("Data/Covariates_dna_methylation.rds")
cpg <- readRDS("Data/Original/imputed_denoised_450K.rds")
proteins <- readRDS("Data/Original/Proteins_selected_denoised_re.rds")
metab <- readRDS(paste0("/Users/barbara/Dropbox/Metabolomics_lung_cancer/Data/Metabolites_imputed_summarised_centroid.rds"))
covars <- readRDS("Data/Original/Covariates_no_rep.rds")

# Changing participant IDs in DNA methylation
covars_cpg <- covars_cpg[which(!is.na(covars_cpg$labnr)), ]
covars_cpg <- covars_cpg[!duplicated(covars_cpg$labnr), ]
cpg <- cpg[which(rownames(cpg) %in% rownames(covars_cpg)), ]
rownames(cpg) <- covars_cpg[rownames(cpg), "labnr"]

# Identifying participants with all measurements
ids <- intersect(rownames(proteins), rownames(metab))
ids <- intersect(ids, rownames(cpg))
cpg <- cpg[ids, ]
proteins <- proteins[ids, ]
metab <- metab[ids, ]
covars <- covars[ids, ]

# Loading the results
stab_proteins <- readRDS("Results/Inflammatory_proteins/Multivariate_lung_cancer.rds")
stab_cpg <- readRDS("Results/DNA_methylation/stability_lung_cancer_all.rds")
stab_metab <- readRDS("Results/Metabolomics/Multivariate_lung_cancer.rds")

# Identifying stable features
ids_proteins <- names(SelectedVariables(stab_proteins))[which(SelectedVariables(stab_proteins) == 1)]
ids_cpg <- names(SelectedVariables(stab_cpg))[which(SelectedVariables(stab_cpg) == 1)]
ids_metab <- names(SelectedVariables(stab_metab))[which(SelectedVariables(stab_metab) == 1)]

# Keeping stable features only
proteins_lc <- proteins[, ids_proteins]
cpg_lc <- cpg[, ids_cpg]
metab_lc <- metab[, ids_metab]

# Saving filtered datasets
dir.create("Data/Filtered_lung_cancer", showWarnings = FALSE)
saveRDS(proteins_lc, "Data/Filtered_lung_cancer/Proteins.rds")
saveRDS(cpg_lc, "Data/Filtered_lung_cancer/DNA_methylation.rds")
saveRDS(metab_lc, "Data/Filtered_lung_cancer/Metabolomics.rds")
saveRDS(covars, "Data/Filtered_lung_cancer/Covariates.rds")
