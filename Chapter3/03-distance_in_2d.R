rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter3/")

library(sharp)

x <- rbind(c(2, 1), c(3, 4))
sqrt(sum((x[1, ] - x[2, ])^2))
dist(x)

pdf(paste0("Working_figures/Distance_2d.pdf"),
  width = 7, height = 7
)
par(mar = c(5, 5, 1, 1), xaxs = "i", yaxs = "i")
plot(x[, 1], x[, 2],
  pch = 3, xlim = c(0, 5), ylim = c(0, 5), cex = 2,
  las = 1, cex.lab = 1.5,
  xlab = "First attribute",
  ylab = "Second attribute",
  bty = "n", lwd = 3
)
text(x[1, 1], x[1, 2], labels = "X1", pos = 2, cex = 1.5, offset = 1)
text(x[2, 1], x[2, 2], labels = "X2", pos = 4, cex = 1.5, offset = 1)
segments(x0 = x[1, 1], y0 = x[1, 2], x1 = x[2, 1], y1 = x[2, 2], col = "red", lwd = 2)
segments(x0 = x[1, 1], y0 = x[1, 2], x1 = x[2, 1], y1 = x[1, 2], col = "blue", lwd = 2)
segments(x0 = x[2, 1], y0 = x[1, 2], x1 = x[2, 1], y1 = x[2, 2], col = "blue", lwd = 2)
points(x[, 1], x[, 2], pch = 3, cex = 2, lwd = 3)
text(mean(x[, 1]), x[1, 2], pos = 1, labels = "a", col = "blue", cex = 1.5)
text(x[2, 1], mean(x[, 2]), pos = 4, labels = "b", col = "blue", cex = 1.5)
text(mean(x[, 1]), mean(x[, 2]), labels = expression(d[12]^E * "=" * sqrt(a^2 + b^2)), col = "red", cex = 1.5, pos = 2)
dev.off()


dist(t(scale(t(x))))
1 - cor(t(x))[1, 2]

simul <- SimulateGraphical()
x <- simul$data

cor(t(x))[1, 2]
1 - (as.matrix(dist(t(scale(t(x)))))[1, 2])^2 / (2 * (ncol(x) - 1))
