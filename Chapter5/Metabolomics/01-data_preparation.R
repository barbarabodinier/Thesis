rm(list = ls())

library(openxlsx)
library(lme4)
library(colorspace)
library(imputeLCMD)
library(DescTools)
source("Scripts/functions.R")

# Creating folders
filepath <- "Figures/1-data_preparation/"
dir.create("Figures", showWarnings = FALSE)
dir.create(filepath, showWarnings = FALSE)
tablepath <- "Tables/1-data_preparation/"
dir.create("Tables", showWarnings = FALSE)
dir.create(tablepath, showWarnings = FALSE)

# Loading the data
covars <- readRDS("Data/Covariates_no_rep.rds")

# Re-calculating packyears from intensity and duration
covars$packyears <- (covars$smok_intensity / 20) * covars$smok_duration
covars$packyears[is.na(covars$packyears)] <- 0
covars$packyears[is.na(covars$smok_intensity)] <- NA
covars$packyears[is.na(covars$smok_duration)] <- NA
range(covars$packyears[which(covars$smoking_status == "Never")])
range(covars$packyears[which(covars$smoking_status == "Former")], na.rm = TRUE)
range(covars$packyears[which(covars$smoking_status == "Current")], na.rm = TRUE)
saveRDS(covars, "Data/Covariates_no_rep.rds")

metab_neg <- read.csv("Data/Original/Lung_cancer_feature_table_RP_NEG_130919.csv")
rownames(metab_neg) <- metab_neg[, 1]
metab_neg <- metab_neg[, -1]
metab_neg <- t(metab_neg)
rownames(metab_neg) <- gsub("^X", "", rownames(metab_neg))
print(dim(metab_neg))

metab_pos <- read.csv("Data/Original/Lung_cancer_feature_table_RP_POS_090819.csv")
rownames(metab_pos) <- metab_pos[, 1]
metab_pos <- metab_pos[, -1]
metab_pos <- t(metab_pos)
rownames(metab_pos) <- gsub("^X", "", rownames(metab_pos))
print(dim(metab_pos))

sample_sheet_neg <- readRDS("Data/Original/Sample_sheet_metabolites_negative_mode.rds")
sample_sheet_pos <- readRDS("Data/Original/Sample_sheet_metabolites_positive_mode.rds")

# Log2-transformation of the data
metab_neg <- log2(metab_neg)
metab_pos <- log2(metab_pos)
table(is.infinite(metab_neg))
table(is.infinite(metab_pos))
metab_neg[is.infinite(metab_neg)] <- NA

# Extracting the blanks
metab_neg_blanks <- metab_neg[grepl("Blank", rownames(metab_neg)), ]
metab_pos_blanks <- metab_pos[grepl("Blank", rownames(metab_pos)), ]

# Extracting the plate and order
for (mode in c("pos", "neg")) {
  sample_sheet <- eval(parse(text = paste0("sample_sheet_", mode)))
  order <- gsub("_.*", "", rownames(sample_sheet))
  names(order) <- rownames(sample_sheet)
  order <- sapply(split(order, f = sample_sheet$Plate), FUN = function(x) {
    mynames <- names(x)
    x <- as.numeric(x)
    x <- x - min(x) + 1
    names(x) <- mynames
    return(x)
  })
  order <- unlist(order)
  names(order) <- gsub(".*\\.", "", names(order))
  plate <- sample_sheet$Plate
  names(plate) <- rownames(sample_sheet)

  assign(paste0("order_all_", mode), order)
  assign(paste0("plate_all_", mode), plate)
}

# Linear mixed models: reference feature
feature_name <- colnames(metab_pos)[grep("474.3188@4.261", colnames(metab_pos))]
metab_pos_samples <- metab_pos[!grepl("Blank", rownames(metab_pos)), ]
x <- metab_pos_samples[, feature_name]
order <- order_all_pos[rownames(metab_pos_samples)]
plate <- plate_all_pos[rownames(metab_pos_samples)]
case <- covars[gsub(".*_", "", names(x)), "case"]
center <- covars[gsub(".*_", "", names(x)), "centre"]
x_blanks <- metab_pos_blanks[, feature_name]

{
  pdf(paste0(filepath, "/Intensity_ref_positive.pdf"), width = 15, height = 7)
  intensityPlot(x = x, x_blanks = x_blanks, plate = plate, case = case, ylim = c(15.5, 16.45))
  dev.off()
}

# Removing the N=7 samples after drop in intensity (positive mode only)
toremove <- sample_sheet_pos$Sample.Name[which(sample_sheet_pos$Sample.Name == 2467243):(which(sample_sheet_pos$Sample.Name == 2467243) + 6)]
print(sample_sheet_pos$Plate[which(sample_sheet_pos$Sample.Name %in% toremove)])
toremove <- rownames(sample_sheet_pos)[which(sample_sheet_pos$Sample.Name %in% toremove)]
metab_pos <- metab_pos[!rownames(metab_pos) %in% toremove, ]
sample_sheet_pos <- sample_sheet_pos[!rownames(sample_sheet_pos) %in% toremove, ]
metab_neg <- metab_neg[!rownames(metab_neg) %in% toremove, ]
sample_sheet_neg <- sample_sheet_neg[!rownames(sample_sheet_neg) %in% toremove, ]

# Extracting the plate and order
for (mode in c("pos", "neg")) {
  sample_sheet <- eval(parse(text = paste0("sample_sheet_", mode)))
  order <- gsub("_.*", "", rownames(sample_sheet))
  names(order) <- rownames(sample_sheet)
  order <- sapply(split(order, f = sample_sheet$Plate), FUN = function(x) {
    mynames <- names(x)
    x <- as.numeric(x)
    x <- x - min(x) + 1
    names(x) <- mynames
    return(x)
  })
  order <- unlist(order)
  names(order) <- gsub(".*\\.", "", names(order))
  plate <- sample_sheet$Plate
  names(plate) <- rownames(sample_sheet)

  assign(paste0("order_all_", mode), order)
  assign(paste0("plate_all_", mode), plate)
}

# Removing the blanks
metab_neg <- metab_neg[!grepl("Blank", rownames(metab_neg)), ]
rownames(metab_neg) <- sapply(lapply(strsplit(rownames(metab_neg), "_"), as.numeric), FUN = function(x) {
  paste(x[1:2], collapse = "_")
})
sample_sheet_neg <- sample_sheet_neg[!grepl("Blank", rownames(sample_sheet_neg)), ]
sample_sheet_neg <- sample_sheet_neg[rownames(metab_neg), ]

metab_pos <- metab_pos[!grepl("Blank", rownames(metab_pos)), ]
rownames(metab_pos) <- sapply(lapply(strsplit(rownames(metab_pos), "_"), as.numeric), FUN = function(x) {
  paste(x[1:2], collapse = "_")
})
sample_sheet_pos <- sample_sheet_pos[!grepl("Blank", rownames(sample_sheet_pos)), ]
sample_sheet_pos <- sample_sheet_pos[rownames(metab_pos), ]

print(all(rownames(metab_neg) == rownames(sample_sheet_neg)))
print(all(rownames(metab_pos) == rownames(sample_sheet_pos)))

# Excluding the features with more than 50% missing
for (mode in c("pos", "neg")) {
  print(mode)
  metab <- eval(parse(text = paste0("metab_", mode)))
  metab_blanks <- eval(parse(text = paste0("metab_", mode, "_blanks")))

  prop_missing <- apply(metab, 2, FUN = function(x) {
    sum(is.na(x))
  }) / nrow(metab)
  print(table(prop_missing < 0.5))
  metab <- metab[, prop_missing < 0.5]
  metab_blanks <- metab_blanks[, prop_missing < 0.5]

  print(ncol(metab))
  assign(paste0("metab_", mode), metab)
  assign(paste0("metab_", mode, "_blanks"), metab_blanks)
}

# Linear mixed models: reference feature
feature_name <- colnames(metab_pos)[grep("474.3188@4.261", colnames(metab_pos))]
x <- metab_pos[, feature_name]
order <- order_all_pos[rownames(metab_pos)]
plate <- plate_all_pos[rownames(metab_pos)]
case <- covars[gsub(".*_", "", names(x)), "case"]
center <- paste0(
  covars[gsub(".*_", "", names(x)), "centre"],
  "_",
  covars[gsub(".*_", "", names(x)), "gender"]
)
id <- gsub(".*_", "", names(x))
mymodel <- lmer(x ~ (order | plate) + (1 | center) + (1 | id))
ranef(mymodel)$plate
ranef(mymodel)$design
x_blanks <- metab_pos_blanks[, feature_name]

# Visualisation
{
  pdf(paste0(filepath, "/Intensity_ref_positive_no_drop_outliers.pdf"), width = 15, height = 7)
  intensityPlot(x = x, x_blanks = x_blanks, plate = plate, case = case, ylim = c(15.5, 16.45))
  dev.off()
}

{
  pdf(paste0(filepath, "/Intensity_ref_positive_no_drop_outliers_slopes.pdf"), width = 15, height = 7)
  intensityPlot(x = x, x_blanks = x_blanks, plate = plate, case = case, model = mymodel, slope = TRUE, ylim = c(15.5, 16.45))
  dev.off()
}

# Denoising of the measurements on plate|order and center
denoised_x <- x - ranef(mymodel)$plate[plate, 1] - ranef(mymodel)$center[center, 1] - ranef(mymodel)$plate[plate, 2] * order # done manually to keep missing data, otherwise: denoised_x=residuals(mymodel)+fixef(mymodel)

# Application of the model to the blanks: denoising only on plate|order (no center)
fixef(mymodel)
mean(metab_pos_blanks[, feature_name], na.rm = TRUE)
denoised_blanks <- metab_pos_blanks[, feature_name] - ranef(mymodel)$plate[plate_all_pos[rownames(metab_pos_blanks)], 1] - ranef(mymodel)$plate[plate_all_pos[rownames(metab_pos_blanks)], 2] * order_all_pos[rownames(metab_pos_blanks)]

# Visualisation
{
  pdf(paste0(filepath, "/Intensity_ref_positive_corrected.pdf"), width = 15, height = 7)
  intensityPlot(x = denoised_x, x_blanks = denoised_blanks, plate = plate, case = case, ylim = c(15.5, 16.45))
  dev.off()
}

# Linear mixed models adjusting on plate, order, centre, gender and id for all features
for (mode in c("pos", "neg")) {
  print(mode)
  metab <- eval(parse(text = paste0("metab_", mode)))
  metab_blanks <- eval(parse(text = paste0("metab_", mode, "_blanks")))
  order_all <- eval(parse(text = paste0("order_all_", mode)))
  plate_all <- eval(parse(text = paste0("plate_all_", mode)))

  denoised <- matrix(NA, ncol = ncol(metab), nrow = nrow(metab))
  rownames(denoised) <- rownames(metab)
  colnames(denoised) <- colnames(metab)
  denoised_blanks <- matrix(NA, ncol = ncol(metab_blanks), nrow = nrow(metab_blanks))
  rownames(denoised_blanks) <- rownames(metab_blanks)
  colnames(denoised_blanks) <- colnames(metab_blanks)

  order <- order_all[rownames(metab)]
  plate <- plate_all[rownames(metab)]
  case <- covars[gsub(".*_", "", rownames(metab)), "case"]
  center <- paste0(
    covars[gsub(".*_", "", rownames(metab)), "centre"],
    "_",
    covars[gsub(".*_", "", rownames(metab)), "gender"]
  )
  id <- gsub(".*_", "", rownames(metab))

  for (feature_name in colnames(metab)) {
    print(which(colnames(metab) == feature_name))
    x <- metab[, feature_name]
    if (sum(duplicated(id[!is.na(x)])) > 0) {
      mymodel <- lmer(x ~ (order | plate) + (1 | center) + (1 | id))
    } else {
      mymodel <- lmer(x ~ (order | plate) + (1 | center))
    }
    x_blanks <- metab_blanks[, feature_name]

    # Denoising of the measurements on plate|order and center
    denoised_x <- x - fixef(mymodel) - ranef(mymodel)$plate[plate, 1] - ranef(mymodel)$center[center, 1] - ranef(mymodel)$plate[plate, 2] * order # done manually to keep missing data
    denoised[, feature_name] <- denoised_x

    # Application of the model to the blanks: denoising only on plate|order (no center)
    denoised_x_blanks <- metab_blanks[, feature_name] - ranef(mymodel)$plate[plate_all[rownames(metab_blanks)], 1] - ranef(mymodel)$plate[plate_all[rownames(metab_blanks)], 2] * order_all[rownames(metab_blanks)]
    denoised_blanks[, feature_name] <- denoised_x_blanks
  }

  assign(paste0("metab_residuals_", mode), denoised)
}

# Identifying outliers
for (mode in c("pos", "neg")) {
  print(mode)
  metab <- eval(parse(text = paste0("metab_residuals_", mode)))

  # Calculating variance over the metabolomics features
  xseq <- apply(metab, 1, var, na.rm = TRUE)

  # Calculating p-value from ANOVA on variance
  thr_list <- sort(unique(xseq), decreasing = TRUE)
  pvalues <- rep(NA, length(thr_list))
  for (k in 1:length(thr_list)) {
    thr <- thr_list[k]
    grouping <- ifelse(xseq >= thr, yes = 1, no = 0)
    myanova <- aov(xseq ~ grouping)
    pvalues[k] <- summary(myanova)[[1]][1, 5]
  }
  grouping <- ifelse(xseq >= thr_list[which.min(pvalues)], yes = 1, no = 0)

  # Storing outputs
  assign(paste0("xseq_", mode), xseq)
  assign(paste0("thr_", mode), thr_list[which.min(pvalues)])
  assign(paste0("outliers_", mode), grouping)
}
detected <- sort(unique(c(
  names(outliers_pos)[which(outliers_pos == 1)],
  names(outliers_neg)[which(outliers_neg == 1)]
)))

{
  pdf(paste0(filepath, "/Detection_outliers.pdf"), width = 15, height = 10)
  par(mar = c(5, 5, 1, 1), mfrow = c(2, 1))
  for (mode in c("pos", "neg")) {
    print(mode)
    xseq <- eval(parse(text = paste0("xseq_", mode)))
    thr <- eval(parse(text = paste0("thr_", mode)))

    mycolours <- rep("navy", length(xseq))
    names(mycolours) <- names(xseq)
    mycolours[detected] <- "red"
    ids <- which(mycolours == "red")
    plot(as.numeric(gsub("_.*", "", names(xseq))), xseq,
      col = mycolours,
      pch = 19, las = 1,
      xlab = "Processing order of the samples",
      ylab = "Residual variance",
      cex.lab = 1.5, cex = 0.7,
      panel.first = c(
        abline(h = 0, lty = 2, col = "black"),
        abline(h = thr, lty = 2, col = "darkred"),
        abline(v = as.numeric(gsub("_.*", "", names(xseq)))[ids], lty = 3, col = lighten("red", amount = 0.5))
      )
    )
    text(as.numeric(gsub("_.*", "", names(xseq)))[ids], xseq[ids],
      labels = 1:length(ids),
      pos = 4, col = "darkred", cex = 0.8, offset = 0.3
    )
  }
  dev.off()
}

# Saving characteristics of outliers
write.table(covars[
  sample_sheet_pos[detected, "Sample.Name"],
  c("age.sample", "gender", "cohort", "centre", "case", "smoking_status")
],
paste0(tablepath, "/Outliers_characteristics.txt"),
row.names = TRUE, quote = FALSE
)
print("Outliers:")
print(length(detected))
print(table((covars[
  sample_sheet_pos[detected, "Sample.Name"],
  "case"
])))

# Excluding outliers
for (mode in c("pos", "neg")) {
  metab <- eval(parse(text = paste0("metab_", mode)))
  metab <- metab[!rownames(metab) %in% detected, ]
  assign(paste0("metab_", mode), metab)
}

# Linear mixed models for denoising of all features
for (mode in c("pos", "neg")) {
  print(mode)
  metab <- eval(parse(text = paste0("metab_", mode)))
  metab_blanks <- eval(parse(text = paste0("metab_", mode, "_blanks")))
  order_all <- eval(parse(text = paste0("order_all_", mode)))
  plate_all <- eval(parse(text = paste0("plate_all_", mode)))

  denoised <- matrix(NA, ncol = ncol(metab), nrow = nrow(metab))
  rownames(denoised) <- rownames(metab)
  colnames(denoised) <- colnames(metab)
  denoised_blanks <- matrix(NA, ncol = ncol(metab_blanks), nrow = nrow(metab_blanks))
  rownames(denoised_blanks) <- rownames(metab_blanks)
  colnames(denoised_blanks) <- colnames(metab_blanks)

  order <- order_all[rownames(metab)]
  plate <- plate_all[rownames(metab)]
  case <- covars[gsub(".*_", "", rownames(metab)), "case"]
  center <- paste0(
    covars[gsub(".*_", "", rownames(metab)), "centre"],
    "_",
    covars[gsub(".*_", "", rownames(metab)), "gender"]
  )
  id <- gsub(".*_", "", rownames(metab))

  for (feature_name in colnames(metab)) {
    print(which(colnames(metab) == feature_name))
    x <- metab[, feature_name]
    # mymodel=lmer(x~(order|plate)+(1|center))
    if (sum(duplicated(id[!is.na(x)])) > 0) {
      mymodel <- lmer(x ~ (order | plate) + (1 | center) + (1 | id))
    } else {
      mymodel <- lmer(x ~ (order | plate) + (1 | center))
    }
    x_blanks <- metab_blanks[, feature_name]

    # Denoising of the measurements on plate|order and center
    denoised_x <- x - ranef(mymodel)$plate[plate, 1] - ranef(mymodel)$center[center, 1] - ranef(mymodel)$plate[plate, 2] * order # done manually to keep missing data, otherwise: denoised_x=residuals(mymodel)+fixef(mymodel)+id
    denoised[, feature_name] <- denoised_x

    # Application of the model to the blanks: denoising only on plate|order (no center)
    denoised_x_blanks <- metab_blanks[, feature_name] - ranef(mymodel)$plate[plate_all[rownames(metab_blanks)], 1] - ranef(mymodel)$plate[plate_all[rownames(metab_blanks)], 2] * order_all[rownames(metab_blanks)]
    denoised_blanks[, feature_name] <- denoised_x_blanks
  }

  assign(paste0("metab_", mode), denoised)
  assign(paste0("metab_", mode, "_blanks"), denoised_blanks)
}

# Correlation between repeated measurements (negative mode)
dup_ids <- sample_sheet_neg$Sample.Name[duplicated(sample_sheet_neg$Sample.Name)]
{
  pdf(paste0(filepath, "Correlation_repeats_preprocessed_metab_neg_log2.pdf"), width = 14, height = 8)
  par(mar = c(5, 5, 3, 1), mfrow = c(1, 2))
  for (k in 1:2) {
    ids <- rownames(sample_sheet_neg)[sample_sheet_neg$Sample.Name == dup_ids[k]]
    plot(metab_neg[ids[1], ], metab_neg[ids[2], ],
      pch = 19, cex = 0.5,
      xlab = "Repeat 1", ylab = "Repeat 2", cex.lab = 1.5, col = "navy",
      panel.first = abline(0, 1, lty = 3),
      main = paste("Participant", k), cex.main = 2
    )
    # mycor <- cor(metab_neg[ids[1], ], metab_neg[ids[2], ], use = "complete.obs")
    mycor <- as.numeric(CCC(metab_neg[ids[1], ], metab_neg[ids[2], ], na.rm = TRUE)$rho.c[1])
    legend("topleft", legend = paste0("r = ", formatC(mycor, format = "f", digits = 2)), cex = 2, bty = "n")
  }
  dev.off()
}

# Correlation between repeated measurements (positive mode)
dup_ids <- sample_sheet_pos$Sample.Name[duplicated(sample_sheet_pos$Sample.Name)]
{
  pdf(paste0(filepath, "Correlation_repeats_preprocessed_metab_pos_log2.pdf"), width = 14, height = 8)
  par(mar = c(5, 5, 3, 1), mfrow = c(1, 2))
  for (k in 1:2) {
    ids <- rownames(sample_sheet_pos)[sample_sheet_pos$Sample.Name == dup_ids[k]]
    plot(metab_pos[ids[1], ], metab_pos[ids[2], ],
      pch = 19, cex = 0.5,
      xlab = "Repeat 1", ylab = "Repeat 2", cex.lab = 1.5, col = "red",
      panel.first = abline(0, 1, lty = 3),
      main = paste("Participant", k), cex.main = 2
    )
    # mycor <- cor(metab_pos[ids[1], ], metab_pos[ids[2], ], use = "complete.obs")
    mycor <- as.numeric(CCC(metab_pos[ids[1], ], metab_pos[ids[2], ], na.rm = TRUE)$rho.c[1])
    legend("topleft", legend = paste0("r = ", formatC(mycor, format = "f", digits = 2)), cex = 2, bty = "n")
    print(cor(metab_pos[ids[1], ], metab_pos[ids[2], ], use = "complete.obs", method = "spearman"))
  }
  dev.off()
}

# Take the average over the repeated measurements
for (mode in c("pos", "neg")) {
  mydata <- eval(parse(text = paste0("metab_", mode)))
  if (mode == "neg") {
    sample_sheet <- sample_sheet_neg
  } else {
    sample_sheet <- sample_sheet_pos
  }
  mydata_no_dup <- mydata[!rownames(mydata) %in% rownames(sample_sheet)[sample_sheet$Sample.Name %in% dup_ids], ]
  rownames(mydata_no_dup) <- sample_sheet[rownames(mydata_no_dup), "Sample.Name"]
  for (k in 1:2) {
    ids <- rownames(sample_sheet)[sample_sheet$Sample.Name == dup_ids[k]]
    mydata_no_dup <- rbind(mydata_no_dup, apply(mydata[ids, ], 2, mean, na.rm = TRUE))
  }
  rownames(mydata_no_dup)[c(nrow(mydata_no_dup) - 1, nrow(mydata_no_dup))] <- dup_ids
  assign(paste0("metab_", mode), mydata_no_dup)
}

# Background subtraction: t-test
for (mode in c("pos", "neg")) {
  metab <- eval(parse(text = paste0("metab_", mode)))
  metab_blanks <- eval(parse(text = paste0("metab_", mode, "_blanks")))

  mytests <- NULL
  for (k in 1:ncol(metab)) {
    if (sum(is.na(metab_blanks[, k])) == nrow(metab_blanks)) {
      mytests <- c(mytests, 0)
    } else {
      if (sum(!is.na(metab_blanks[, k])) == 1) {
        mytest <- t.test(x = metab[, k], mu = mean(metab_blanks[, k], na.rm = TRUE), alternative = "greater")
        mytests <- c(mytests, mytest$p.value)
      } else {
        mytest <- t.test(x = metab[, k], y = metab_blanks[, k], alternative = "greater")
        mytests <- c(mytests, mytest$p.value)
      }
    }
  }
  print(table(mytests < 0.05))
  metab <- metab[, mytests < 0.05]
  print(ncol(metab))
  assign(paste0("metab_", mode), metab)
  assign(paste0("metab_", mode, "_blanks"), metab_blanks)
}

# Saving the datasets
saveRDS(metab_neg, "Data/Metabolites_negative_no_dup.rds")
saveRDS(metab_pos, "Data/Metabolites_positive_no_dup.rds")

# Saving imputed datasets
for (mode in c("pos", "neg")) {
  print(mode)
  metab <- eval(parse(text = paste0("metab_", mode)))

  # Data imputation
  imputed <- impute.QRILC(metab)[[1]]

  # Exclusion of participants with missing bmi or packyears
  print("Missing covariates:")
  covars <- readRDS("Data/Covariates_no_rep.rds")
  covars <- covars[rownames(metab), ]
  ids <- which(is.na(covars$bmi) | is.na(covars$packyears))
  print(length(ids))
  print(table(covars[ids, "case"]))
  imputed <- imputed[which(!is.na(covars$bmi) & !is.na(covars$packyears)), ]

  if (mode == "pos") {
    saveRDS(imputed, paste0("Data/Metabolites_positive_imputed.rds"))
  } else {
    saveRDS(imputed, paste0("Data/Metabolites_negative_imputed.rds"))
  }
}
