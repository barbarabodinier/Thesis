rm(list = ls())

library(openxlsx)
library(colorspace)
library(basicPlotteR)
library(plotrix)
library(sharp)

# Parameters
tmpcol <- "navy"

for (subtype in c("lung_cancer", "small-cell_carcinoma", "adenocarcinoma")) {
  print(subtype)

  # Loading the results
  tmpstab <- readRDS(paste0("Results/Inflammatory_proteins/Multivariate_", subtype, ".rds"))
  tmpexpl <- readRDS(paste0("Results/Inflammatory_proteins/Incremental_performances_", subtype, ".rds"))
  tmpsmoking <- readRDS("Results/Inflammatory_proteins/Univariate_packyears.rds")

  # Calibration and selection proportions
  pdf(paste0("Working_figures/Inflammatory_proteins/Stability_selection_", subtype, ".pdf"),
    width = 20, height = 6
  )
  layout(mat = rbind(c(1, 2)), widths = c(1, 2))
  par(mar = c(7, 5, 5, 7))
  CalibrationPlot(tmpstab)
  selprop <- sort(SelectionProportions(tmpstab), decreasing = TRUE)
  plot(selprop,
    type = "h", lend = 1, lwd = 5,
    xlab = "", xaxt = "n",
    col = ifelse(SelectedVariables(tmpstab)[names(selprop)] == 1,
      yes = tmpcol, no = lighten(tmpcol, amount = 0.85)
    ),
    ylab = "Selection proportion", cex.lab = 1.5,
    ylim = c(0, 1),
    panel.first = abline(h = Argmax(tmpstab)[2], lty = 2, col = "darkred")
  )
  for (k in 1:length(selprop)) {
    axis(
      side = 1, at = k, labels = names(selprop)[k],
      las = 2, cex.axis = 0.8,
      col.axis = darken(ifelse(SelectedVariables(tmpstab)[names(selprop)] == 1,
        yes = tmpcol, no = lighten(tmpcol, amount = 0.85)
      )[k], amount = 0.2)
    )
  }
  dev.off()

  # Classifying into smoking mediators and independent
  selprop <- SelectionProportions(tmpstab)
  pvalues <- -log10(tmpsmoking[, 3]) * sign(tmpsmoking[, 2])
  bonf_thr <- 0.05 / nrow(tmpsmoking)
  pi_thr <- Argmax(tmpstab)[2]
  ids_signif <- rownames(tmpsmoking)[p.adjust(tmpsmoking[, 3], method = "bonf") < 0.05]
  ids_selected <- names(SelectedVariables(tmpstab))[which(SelectedVariables(tmpstab) == 1)]
  ids_highlight <- unique(c(ids_signif, ids_selected))
  names(pvalues) <- rownames(tmpsmoking)
  selprop <- selprop[names(pvalues)]

  pdf(paste0("Working_figures/Inflammatory_proteins/Smoking_and_disease_", subtype, ".pdf"),
    width = 7, height = 7
  )
  par(mfrow = c(1, 1), mar = c(5, 5, 1, 1))
  plot(NA,
    xlim = range(pvalues), ylim = c(0, 1),
    xlab = expression(-log[10] * "(p-value of packyears) (signed)"),
    ylab = paste0("Selection proportion for ", gsub("_", " ", subtype)),
    cex.lab = 1.5, las = 1,
    panel.first = c(
      abline(v = 0, lty = 3),
      abline(h = Argmax(tmpstab)[2], lty = 2, col = "darkred"),
      abline(v = -log10(bonf_thr), lty = 2, col = "tan"),
      abline(v = log10(bonf_thr), lty = 2, col = "tan")
    )
  )
  polygon(
    x = c(min(pvalues) * 2, log10(bonf_thr), log10(bonf_thr), min(pvalues) * 2),
    y = c(-1, -1, 2, 2), border = NA,
    col = adjustcolor(col = "tan", alpha.f = 0.1)
  )
  polygon(
    x = c(max(pvalues) * 2, -log10(bonf_thr), -log10(bonf_thr), max(pvalues) * 2),
    y = c(-1, -1, 2, 2), border = NA,
    col = adjustcolor(col = "tan", alpha.f = 0.1)
  )
  polygon(
    x = c(min(pvalues) * 2, max(pvalues) * 2, max(pvalues) * 2, min(pvalues) * 2),
    y = c(pi_thr, pi_thr, 2, 2), border = NA,
    col = adjustcolor(col = "red", alpha.f = 0.1)
  )
  points(pvalues, selprop,
    pch = 19,
    col = lighten(tmpcol, amount = 0.85),
    panel.first = c(
      abline(v = 0, lty = 3),
      abline(h = Argmax(tmpstab)[2], lty = 2, col = "darkred"),
      abline(v = -log10(0.05 / nrow(tmpsmoking)), lty = 2, col = "tan"),
      abline(v = log10(0.05 / nrow(tmpsmoking)), lty = 2, col = "tan")
    )
  )
  points(pvalues[ids_selected],
    selprop[ids_selected],
    pch = 19, col = tmpcol
  )
  addTextLabels(pvalues[ids_highlight],
    selprop[ids_highlight],
    keepLabelsInside = TRUE,
    col.label = ifelse(ids_highlight %in% ids_selected,
      yes = darken(tmpcol, amount = 0.2),
      no = lighten(tmpcol, amount = 0.6)
    ),
    col.line = ifelse(ids_highlight %in% ids_selected,
      yes = tmpcol,
      no = lighten(tmpcol, amount = 0.85)
    ),
    col.background = NA,
    border = NA,
    cex.label = 1,
    labels = ids_highlight
  )
  dev.off()

  # Explanatory performances
  plotname <- paste0("Working_figures/Inflammatory_proteins/Explanatory_performances_", subtype, ".pdf")
  pdf(plotname, width = 20, height = 6)
  layout(mat = rbind(c(1, 2)), widths = c(1, 2))
  par(mar = c(10, 5, 5, 7))
  plot.new()
  tmpexpl$names[1] <- "Packyears"
  IncrementalPlot(tmpexpl,
    # ylim = c(0.65, 0.85),
    ylab = "AUC", ylas = 0, sfrac = 0.002,
    col = c(
      "darkred",
      rep(tmpcol, sum(SelectedVariables(tmpstab))),
      rep(
        darken(lighten(tmpcol, amount = 0.85), amount = 0.2),
        ncol(tmpstab$params$xdata) - sum(SelectedVariables(tmpstab))
      )
    ),
    quantiles = c(0.25, 0.75), xcex.axis = 0.8, xgrid = TRUE
  )
  abline(h = median(tmpexpl$AUC[[sum(SelectedVariables(tmpstab) == 1) + 1]]), lty = 2, col = tmpcol)
  for (k in 2:length(tmpexpl$names)) {
    axis(
      side = 3, at = k, las = 2, cex.axis = 0.8,
      labels = format(selprop[tmpexpl$names], format = "f", digits = 2)[k],
      col.axis = c(
        "darkred",
        rep(tmpcol, sum(SelectedVariables(tmpstab))),
        rep(
          darken(lighten(tmpcol, amount = 0.85), amount = 0.2),
          ncol(tmpstab$params$xdata) - sum(SelectedVariables(tmpstab))
        )
      )[k]
    )
  }
  dev.off()
  system(paste("pdfcrop --margin 10", plotname, plotname))
}
