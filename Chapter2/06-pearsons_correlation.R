rm(list = ls())

library(MASS)

# Simulation parameters
n <- 1000
p <- 2

{
  pdf("Working_figures/Pearsons_correlations.pdf",
    width = 25, height = 5
  )
  par(mfrow = c(1, 5), mar = c(5, 5, 4, 2))
  for (sigma_xy in c(-0.9, -0.5, 0, 0.5, 0.9)) {

    # Simulating data from multivariate normal distribution
    set.seed(1)
    sigma <- matrix(c(1, sigma_xy, sigma_xy, 1), ncol = 2, byrow = TRUE)
    x <- MASS::mvrnorm(n, rep(0, p), sigma)
    colnames(x) <- paste0("var", 1:ncol(x))
    rownames(x) <- paste0("obs", 1:nrow(x))

    # Making figure
    plot(x[, 1], x[, 2],
      las = 1, pch = 19, col = "navy",
      xlab = expression(X[1]), ylab = expression(X[2]),
      cex.lab = 2,
      panel.first = c(
        abline(h = 0, lty = 2),
        abline(v = 0, lty = 2)
      )
    )
    mtext(
      text = eval(parse(text = paste0(
        "expression(rho*'=",
        formatC(sigma_xy, format = "f", digits = 3),
        " ('*hat(rho)*'=",
        formatC(cor(x[, 1], x[, 2]), format = "f", digits = 3), ")')"
      ))),
      side = 3, cex = 2, line = 0.25
    )
  }
  dev.off()
}
