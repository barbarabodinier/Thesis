rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Latest/chapter4/")

library(fake)
library(corpcor)
library(abind)

# Different proportions of explained variances
set.seed(1)
theta <- SimulateAdjacency(pk = 20, nu_within = 0.8)
ev_list <- c(0.2, 0.4, 0.6)

pdf("Working_figures/Correlations_different_ev.pdf",
  width = 15, height = 4.5
)
par(mar = c(4, 4, 1, 6), mfrow = c(1, 3))
for (k in 1:3) {
  set.seed(1)
  omega <- SimulatePrecision(
    theta = theta,
    pd_strategy = "min_eigenvalue",
    ev = ev_list[k]
  )
  Heatmap(
    mat = cov2cor(solve(omega$omega)),
    col = c("navy", "white", "darkred"),
    legend_range = c(-1, 1),
    legend = ifelse(k == 3, yes = TRUE, no = FALSE)
  )
}
dev.off()

# Maximisation of the contrast
set.seed(1)
omega <- SimulatePrecision(
  theta = theta,
  pd_strategy = "min_eigenvalue"
)
omega_tilde <- omega$omega

# Storing contrast for different values of u
u_list <- 10^-(seq(-1, 5, by = 0.1))
contrasts <- rep(NA, length(u_list))
mycor <- NULL
diag(omega_tilde) <- 0
for (k in 1:length(u_list)) {
  u <- u_list[k]

  # Making positive definite with given u
  set.seed(1)
  omega <- MakePositiveDefinite(
    omega = omega_tilde,
    pd_strategy = "min_eigenvalue",
    u_list = u,
    ev_xx = NULL
  )

  # Calculating correlation matrix
  mycor <- abind(mycor, cov2cor(solve(omega$omega)), along = 3)

  # Calculating contrast
  contrasts[k] <- Contrast(cov2cor(solve(omega$omega)))
}

plotname <- "Working_figures/Correlations_contrast.pdf"
pdf(plotname, width = 12, height = 8)
layout(mat = matrix(c(1, 1, 1, 2, 3, 4), ncol = 3, byrow = TRUE), heights = c(1.25, 1))
par(mar = c(5, 5, 1, 5))

# Contrast as a function of u
par(xpd = FALSE)
plot(log(u_list), contrasts,
  pch = 19, col = "navy", cex = 0.5,
  xlab = "log(u)", ylab = "Contrast", cex.lab = 1.5, las = 1, bty = "n"
)
abline(v = range(log(u_list)), lty = 2, col = "black")
abline(v = axTicks(1), lty = 3, col = "grey")
abline(h = axTicks(2), lty = 3, col = "grey")
points(log(u_list), contrasts, pch = 19, col = "navy", cex = 0.5)
points(log(u_list)[which.max(contrasts)], max(contrasts), pch = 19, col = "red")
abline(v = log(u_list)[which.max(contrasts)], lty = 2, col = "red")
abline(h = max(contrasts), lty = 2, col = "red")

# Corresponding correlation matrices
Heatmap(mycor[, , dim(mycor)[3]],
  col = c("navy", "white", "darkred"),
  legend_range = c(-1, 1), legend = FALSE, axes = FALSE
)
Heatmap(mycor[, , which.max(contrasts)],
  col = c("navy", "white", "darkred"),
  legend_range = c(-1, 1), legend = FALSE, axes = FALSE
)
Heatmap(mycor[, , 1],
  col = c("navy", "white", "darkred"),
  legend_range = c(-1, 1), legend_length = 10, legend = TRUE, axes = FALSE
)
dev.off()
system(paste("pdfcrop --margin 10", plotname, plotname))
