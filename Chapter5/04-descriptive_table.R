rm(list = ls())

library(tableone)

# Loading the data
proteins <- readRDS("Data/Original/Proteins_selected_denoised_re.rds")
metab <- readRDS("Data/Original/Metabolites_positive_imputed.rds")
cpg <- readRDS("Data/Original/imputed_denoised_450K.rds")
covars <- readRDS("Data/Original/Covariates_no_rep.rds")
covars_cpg_dp <- readRDS("Data/Original/covariates_cpg.rds")

# Creating the gender_cohort variable
covars$gender_cohort <- paste0(covars$cohort, "_", covars$gender)

# Setting NOWAC as first centre
covars$centre <- factor(covars$centre)
covars$centre <- relevel(covars$centre, ref = "NOWAC")

# Setting other lung cancer as last subtype
subtype_list <- sort(unique(covars$subtype[!is.na(covars$subtype)]))
covars$subtype <- factor(covars$subtype, levels = c(
  "Small-cell carcinoma",
  "Adenocarcinoma",
  "Squamous-cell carcinoma",
  "Large-cell carcinoma",
  "Other lung cancer"
))

# Preparing methylation-specific dataset
covars_cpg <- covars_cpg_dp[rownames(cpg), ]
covars_cpg <- covars_cpg[which(!is.na(covars_cpg$bmi) & !is.na(covars_cpg$packyears)), ]
ids_cpg <- covars_cpg$labnr[which(covars_cpg$lc_bc == "lc")]
ids_sample <- rownames(covars_cpg)[which(covars_cpg$lc_bc == "lc")]
covars_cpg <- rbind(
  covars_cpg[ids_sample[which(ids_cpg %in% rownames(covars))], ],
  covars_cpg[which(covars_cpg$lc_bc == "bc"), ]
)

# Creating the gender_cohort variable
covars_cpg$gender_cohort <- paste0(covars_cpg$cohort, "_", covars_cpg$gender)
levels(covars_cpg$center) <- c("Florence", "Naples", "Ragusa", "Turin", "Varese", "NOWAC")
covars_cpg$centre <- covars_cpg$center

# Setting NOWAC as first centre
covars_cpg$centre <- relevel(covars_cpg$centre, ref = "NOWAC")

# Preparing protein-specific dataset
covars_proteins <- covars[rownames(proteins), ]
covars_proteins <- covars_proteins[which(!is.na(covars_proteins$bmi) & !is.na(covars_proteins$packyears)), ]

# Preparing metabolomics-specific dataset
covars_metab <- covars[rownames(metab), ]

# Preparing integration-specific dataset
covars_cpg_dp <- covars_cpg_dp[which(!is.na(covars_cpg_dp$labnr)), ]
covars_cpg_dp <- covars_cpg_dp[!duplicated(covars_cpg_dp$labnr), ]
cpg <- cpg[which(rownames(cpg) %in% rownames(covars_cpg_dp)), ]
rownames(cpg) <- covars_cpg_dp[rownames(cpg), "labnr"]
ids <- intersect(rownames(proteins), rownames(metab))
ids <- intersect(ids, rownames(cpg))
covars_integration <- covars[ids, ]

# Adding subtypes to CpG data
covars_cpg$subtype <- covars[covars_cpg$labnr, "subtype"]
covars_cpg$subtype[which(covars_cpg$lc_bc == "bc")] <- NA
table(covars_cpg$case_lc, covars_cpg$subtype, useNA = "always")
covars_cpg$subtype[which(is.na(covars_cpg$subtype) & (covars_cpg$case_lc == 1))] <- "Unknown"
covars_cpg$subtype <- factor(covars_cpg$subtype, levels = levels(covars$subtype))

# Setting ttd to missing in controls
covars_cpg$ttd[which(covars_cpg$case_lc == 0)] <- NA

# Checking numbers
print(nrow(covars_cpg))
print(nrow(covars_proteins))
print(nrow(covars_metab))
print(nrow(covars_integration))

# Preparing table one
mytableone <- NULL
for (study_name in c("cpg", "proteins", "metab", "integration")) {
  print(study_name)
  tmpcovars <- eval(parse(text = paste0("covars_", study_name)))
  if (study_name == "cpg") {
    myvars <- c(
      "age.recr", "bmi", "gender_cohort", "centre",
      "smoking_status", "smoking_duration", "smoking_intensity", "packyears",
      "ttd", "subtype"
    )
    mystrata <- "case_lc"
  } else {
    myvars <- c(
      "age.sample", "bmi", "gender_cohort", "centre",
      "smoking_status", "smok_duration", "smok_intensity", "packyears",
      "ttd", "subtype"
    )
    mystrata <- "case"
  }
  out <- suppressWarnings(CreateTableOne(
    vars = myvars,
    data = tmpcovars,
    strata = mystrata,
    addOverall = FALSE
  ))
  continuous_pvalues <- attributes(out$ContTable)$pValues
  categorical_pvalues <- attributes(out$CatTable)$pValues
  tmptableone <- print(out,
    printToggle = FALSE, noSpaces = TRUE,
    catDigits = 2, contDigits = 2, explain = FALSE
  )
  tmptableone[rownames(continuous_pvalues), "p"] <- format(continuous_pvalues$pNormal, format = "e", digits = 2)
  tmptableone[rownames(categorical_pvalues), "p"] <- format(categorical_pvalues$pApprox, format = "e", digits = 2)
  tmptableone <- tmptableone[, 1:3]
  mytableone <- cbind(mytableone, tmptableone)
}

# Cleaning table
mytableone[grep("NA", mytableone)] <- ""
mytableone[grep("NaN", mytableone)] <- ""
mytableone <- gsub("( ", "(", mytableone, fixed = TRUE)
mytableone <- gsub("1.0e+00", ">0.99", mytableone, fixed = TRUE)
rownames(mytableone) <- gsub("^   ", "", rownames(mytableone))

# Renaming variables
rownames(mytableone)[which(rownames(mytableone) == "age.recr")] <- "Age at blood collection (years)"
rownames(mytableone)[which(rownames(mytableone) == "bmi")] <- "Body Mass Index (kg/m2)"
rownames(mytableone)[which(rownames(mytableone) == "gender_cohort")] <- "Sex and cohort"
rownames(mytableone)[which(rownames(mytableone) == "centre")] <- "Recruitment centre"
rownames(mytableone)[which(rownames(mytableone) == "smoking_status")] <- "Smoking status"
rownames(mytableone)[which(rownames(mytableone) == "smoking_duration")] <- "Smoking duration (years)"
rownames(mytableone)[which(rownames(mytableone) == "smoking_intensity")] <- "Smoking intensity (cigarettes/day)"
rownames(mytableone)[which(rownames(mytableone) == "packyears")] <- "Packyears"
rownames(mytableone)[which(rownames(mytableone) == "ttd")] <- "Time-to-diagnosis (years)"
rownames(mytableone)[which(rownames(mytableone) == "subtype")] <- "Lung cancer subtype"

# Saving prepared table
write.table(mytableone, paste0("Tables/Table1.txt"),
  row.names = TRUE, col.names = TRUE,
  quote = FALSE, eol = "££\n", sep = "&"
)

# Saving covariates for DNA methylation
saveRDS(covars_cpg, "Data/Covariates_dna_methylation.rds")
