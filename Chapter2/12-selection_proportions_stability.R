rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter2/")

library(glmnet)
library(sharp)

source("Scripts/additional_functions_specific_to_comparisons.R")

# Data simulation
set.seed(0)
simul <- SimulateRegression(
  n = 100, pk = 50, nu_within = 0,
  nu_xz = 0.2, ev_xz = 0.8
)

stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata)
Lambda=stab$Lambda[c(2, 5, ArgmaxId(stab)[1], 35, 50)]
stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata, Lambda=Lambda)

ArgmaxId(stab)

pdf("Working_figures/Selection_proportions_stability.pdf",
    width = 12, height = 14
)
par(mfrow = c(5, 1), mar = c(3, 5, 1, 1))
for (i in 1:5) {
  plot(SelectionProportions(stab, argmax_id = i),
       type = "h", lend = 1, lwd = 10, col = "navy",
       xlab = "",
       ylab = eval(parse(text = paste0(
         "expression(lambda*' = ",
         formatC(stab$Lambda[i], format = "e", digits = 2),
         " (q = ", stab$Q[i], ") ')"
       ))),
       cex.lab = 1.5, xaxt = "n", ylim = c(0, 1),
       panel.first = c(
         abline(h = 0, lty = 2),
         abline(v = seq(1, ncol(simul$xdata)), lty = 3, col = "grey")
       ),
       las = 1
  )
  axis(side=2, at=1, labels=LETTERS[i], cex.axis=4, line=1.5, tick=FALSE, las=2)
}
for (k in 1:ncol(simul$xdata)) {
  axis(side = 1, at = k, labels = paste0("X", k), las = 2)
}
dev.off()

pdf("Working_figures/Calibration_plot_example_stability.pdf",
    width = 7, height = 10
)
par(mar=rep(7, 4))
CalibrationPlot(stab)
dev.off()

