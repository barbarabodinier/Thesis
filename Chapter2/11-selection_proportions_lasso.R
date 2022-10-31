rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter2/")

library(glmnet)
library(sharp)

source("Scripts/additional_functions_specific_to_comparisons.R")

# Data simulation
set.seed(1)
simul <- SimulateRegression(
  n = 100, pk = 50, nu_within = 0,
  nu_xz = 0.2, ev_xz = 0.5
)

set.seed(1)
out <- cv.glmnet(x = simul$xdata, y = simul$ydata)
q <- sum(coef(out, s = "lambda.1se")[-1] != 0)

stab <- VariableSelection(xdata = simul$xdata, ydata = simul$ydata)
CalibrationPlot(stab)

id <- which.min(abs(stab$Q - q))
stab$Lambda[id]

beta1 <- stab$Beta[id, , 1]
beta2 <- stab$Beta[id, , 3]

pdf("Working_figures/Beta_coefficients_lasso.pdf",
  width = 12, height = 9
)
par(mar = c(5, 5, 1, 1), mfrow = c(3, 1))
plot(beta1,
  type = "h", lend = 1, lwd = 10, col = "darkolivegreen",
  xlab = "", ylab = expression(beta * "-coefficient (Run 1)"),
  cex.lab = 1.5, xaxt = "n",
  panel.first = c(
    abline(h = 0, lty = 2),
    abline(v = seq(1, length(beta1)), lty = 3, col = "grey")
  ),
  las = 1
)
for (k in 1:length(beta1)) {
  axis(
    side = 1, at = k, labels = paste0("X", k), las = 2,
    col.axis = ifelse(beta1[k] != 0, yes = "darkolivegreen", no = "grey")
  )
}
plot(beta2,
  type = "h", lend = 1, lwd = 10, col = "tomato",
  xlab = "", ylab = expression(beta * "-coefficient (Run 2)"),
  cex.lab = 1.5, xaxt = "n",
  panel.first = c(
    abline(h = 0, lty = 2),
    abline(v = seq(1, length(beta1)), lty = 3, col = "grey")
  ),
  las = 1
)
for (k in 1:length(beta1)) {
  axis(
    side = 1, at = k, labels = paste0("X", k), las = 2,
    col.axis = ifelse(beta2[k] != 0, yes = "tomato", no = "grey")
  )
}
plot(SelectionProportions(stab, argmax_id = id),
  type = "h", lend = 1, lwd = 10, col = "navy",
  xlab = "", ylab = "Selection proportion",
  cex.lab = 1.5, xaxt = "n", ylim = c(0, 1),
  panel.first = c(
    abline(h = 0, lty = 2),
    abline(v = seq(1, length(beta1)), lty = 3, col = "grey")
  ),
  las = 1
)
for (k in 1:length(beta2)) {
  axis(side = 1, at = k, labels = paste0("X", k), las = 2)
}
dev.off()
