rm(list = ls())

pdf("Working_figures/Normal_distribution.pdf")
xseq <- seq(-5, 5, length.out = 1e5)
par(mar = c(5, 5, 1, 1))
plot(xseq, dnorm(x = xseq),
  type = "l", col = "navy", las = 1,
  xlab = expression(x), ylab = "Density", cex.lab = 1.5
)
abline(v = 0, lty = 2, col = "red")
abline(v = 1, lty = 3, col = "forestgreen")
abline(v = -1, lty = 3, col = "forestgreen")
dev.off()
