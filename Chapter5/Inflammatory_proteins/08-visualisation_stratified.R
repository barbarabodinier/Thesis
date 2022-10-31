rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter5/")

library(openxlsx)
library(colorspace)
library(basicPlotteR)
library(plotrix)
library(sharp)

# Defining strata
stratum_name <- c("Never smokers", "Short time-to-diagnosis", "Long time-to-diagnosis", "Men")
names(stratum_name) <- c("never_smokers", "small_ttd", "large_ttd", "men")

for (stratum in names(stratum_name)) {
  print(stratum)

  # Loading results on full sample
  if (stratum == "never_smokers") {
    subtype <- "adenocarcinoma"
  } else {
    subtype <- "lung_cancer"
  }
  stab <- readRDS(paste0("Results/Inflammatory_proteins/Multivariate_", subtype, ".rds"))
  expl <- readRDS(paste0("Results/Inflammatory_proteins/Incremental_performances_", subtype, ".rds"))
  smoking <- readRDS("Results/Inflammatory_proteins/Univariate_packyears.rds")

  # Loading results from stratified analyses
  tmpstab <- readRDS(paste0("Results/Inflammatory_proteins/Multivariate_lung_cancer_", stratum, ".rds"))
  tmpexpl <- readRDS(paste0("Results/Inflammatory_proteins/Incremental_performances_lung_cancer_", stratum, ".rds"))

  ids_highlight <- sort(unique(c(
    names(SelectedVariables(stab))[which(SelectedVariables(stab) == 1)],
    names(SelectedVariables(tmpstab))[which(SelectedVariables(tmpstab) == 1)]
  )))

  pdf(paste0("Working_figures/Inflammatory_proteins/Scatter_plot_", stratum, ".pdf"),
    width = 14, height = 7
  )
  par(mar = c(7, 5, 5, 7), mfrow = c(1, 2))
  tmpcol <- "navy"
  CalibrationPlot(tmpstab)
  plot(SelectionProportions(stab),
    SelectionProportions(tmpstab),
    xlab = ifelse(stratum == "men", yes = "Women", no = "Full sample"),
    ylab = stratum_name[stratum],
    cex.lab = 1.5,
    pch = 19, las = 1,
    xlim = c(0, 1), ylim = c(0, 1),
    col = lighten(tmpcol, amount = 0.85),
    panel.first = c(
      abline(v = Argmax(stab)[2], col = "darkred", lty = 2),
      abline(h = Argmax(tmpstab)[2], col = "darkred", lty = 2)
    )
  )
  points(SelectionProportions(stab)[ids_highlight],
    SelectionProportions(tmpstab)[ids_highlight],
    pch = 19, col = tmpcol
  )
  addTextLabels(SelectionProportions(stab)[ids_highlight],
    SelectionProportions(tmpstab)[ids_highlight],
    keepLabelsInside = TRUE,
    col.label = darken(tmpcol, amount = 0.2),
    col.line = tmpcol,
    col.background = NA,
    border = NA,
    cex.label = 1,
    labels = ids_highlight
  )
  dev.off()
}
