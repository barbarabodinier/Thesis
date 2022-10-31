rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter4/")

library(fake)
library(sharp)
library(igraph)

source("Scripts/additional_functions_specific_to_comparisons_stability.R")

# Simulation parameters
simul_study_id <- 1
topology <- "random"

# Template design
pi_list <- seq(0.6, 0.9, by = 0.05)
pal_subsampling <- colorRampPalette(c("darkred", "tomato"))(9)
mycolours <- lighten(c(pal_subsampling, "navy", "darkolivegreen"), amount = 0.2)
mypch <- c(rep(17, length(pi_list)), rep(17, length(pi_list)), 18, 18, 18, 18)
dimensionality <- c("Low", "Intermediate", "High")

# Saving figure
for (metric in c("F1_score", "Precision", "Recall")) {
  pdf(paste0("Working_figures/Boxplot_", tolower(metric), "_tau_", simul_study_id, "_", topology, ".pdf"),
    width = 12, height = 4.5
  )
  par(mar = c(7, 5, 3, 1), mfrow = c(1, 3))
  for (simul_id in 1:3) {
    performances <- readRDS(paste0("Results/Stability_selection/2-simulations/1-graphical_model/Sensitivity_tau_", simul_study_id, "_", topology, "/Performances_", simul_id, "_merged.rds"))
    performances <- performances[, , 1:1000]
    colnames(performances) <- tolower(colnames(performances))
    mylist <- list()
    for (k in 1:nrow(performances)) {
      mylist <- c(mylist, list(as.numeric(performances[k, tolower(metric), ])))
    }
    xseq <- c(1:9, 11, 13)
    myylab <- ifelse(simul_id == 1, yes = metric, no = "")
    if (myylab == "F1_score") {
      myylab <- eval(parse(text = "expression(F[1]*'-score')"))
    }
    boxplot(
      at = xseq, mylist, col = mycolours,
      boxcol = "white", whiskcol = mycolours, staplecol = mycolours, medcol = darken(mycolours, amount = 0.5),
      whisklty = 1, range = 0, las = 2, main = dimensionality[simul_id], cex.main = 2,
      ylab = myylab, cex.lab = 1.5, xaxt = "n", ylim = c(0, 1), boxwex = 0.5
    )
    abline(h = axTicks(2), lty = 3, col = "grey")
    boxplot(
      at = xseq, mylist, col = mycolours,
      boxcol = "white", whiskcol = mycolours, staplecol = mycolours, medcol = darken(mycolours, amount = 0.5),
      whisklty = 1, range = 0, las = 2, add = TRUE,
      ylab = myylab, cex.lab = 1.5, xaxt = "n", boxwex = 0.5
    )
    axis(side = 1, at = xseq, labels = c(seq(0.1, 0.9, by = 0.1), "CPSS", "Bootstrapping"), las = 2)
    axis(side = 1, at = c(0.5, 9.5), labels = NA, line = 3)
    axis(side = 1, at = mean(c(0.5, 9.5)), labels = "Subsampling", line = 3, tick = FALSE)
    abline(v = c(10, 12), lty = 3, col = "black")
  }
  dev.off()
}
