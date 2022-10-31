rm(list = ls())

library(fake)
library(sharp)
library(MASS)
library(colorspace)
library(mnormt)
library(cluster)

source("Scripts/additional_functions_specific_to_comparisons.R")

# Simulation of data with clusters
set.seed(0)
n <- c(20, 20, 20)
simul <- SimulateClustering(
  n = n,
  pk = 2,
  theta_xc = c(1, 1),
  ev_xc = 0.6,
  output_matrices = TRUE
)
x <- simul$data

pdf(paste0("Working_figures/Silhouette_width_true_clustering.pdf"),
  width = 8, height = 12
)
mycolours <- c("royalblue", "tomato", "tan")
layout(mat = cbind(c(1, 2)), heights = c(2, 1))
par(mar = c(5, 5, 1, 1))
ids <- c(589, 229, 106)
plot(simul$data,
  pch = 19, cex = 3,
  col = mycolours[simul$theta],
  xlab = "First attribute",
  ylab = "Second attribute",
  cex.lab = 1.5, las = 1
)
text(simul$data, labels = 1:nrow(simul$data))

for (k in 1:nrow(simul$mu_mixture)) {
  x <- seq(-5, 5, length.out = 1000)
  y <- seq(-5, 5, length.out = 1000)
  mu <- simul$mu_mixture[k, ]
  sigma <- simul$sigma
  f <- function(x, y) dmnorm(cbind(x, y), mu, sigma)
  z <- outer(x, y, f)

  contour(x, y, z,
    las = 1, cex.main = 1.5,
    nlevels = 10,
    xlab = "", ylab = "", bty = "n",
    lwd = 2, cex.lab = 1.5,
    col = darken(mycolours[k], amount = 0.5),
    pch = 19, cex = 0.5, las = 1, lty = 1,
    drawlabels = TRUE, add = TRUE
  )
}

points(simul$data,
  pch = 19, cex = 3,
  col = mycolours[simul$theta],
  xlab = "X1", ylab = "X2", cex.lab = 1.5, las = 1
)
text(simul$data, labels = 1:nrow(simul$data))

true_sil <- silhouette(x = simul$theta, dist = dist(simul$data))[, 3]
plot(true_sil,
  type = "h", lwd = 5, lend = 1, bty = "n",
  las = 1, xlab = "", xaxt = "n", ylim = c(-0.4, 0.7),
  ylab = "Silhouette width", cex.lab = 1.5,
  col = mycolours[simul$theta],
  panel.first = abline(v = 1:nrow(simul$data), lty = 2, col = mycolours[simul$theta])
)
for (k in 1:nrow(simul$data)) {
  axis(
    side = 1, at = k, labels = k, las = 2, col.axis = mycolours[simul$theta[k]],
    tick = FALSE, line = -0.5
  )
}
for (k in 1:length(unique(simul$theta))) {
  axis(side = 1, at = range(which(simul$theta == k)) + c(-0.25, 0.25), labels = NA, line = 2)
  axis(
    side = 1, at = mean(which(simul$theta == k)), labels = paste0("Simulated cluster ", k),
    line = 1.5, col.axis = darken(mycolours[k], amount = 0.5), tick = FALSE
  )
}
dev.off()


pdf(paste0("Working_figures/Silhouette_width_pam_clustering.pdf"),
  width = 8, height = 12
)
pam <- pam(x = simul$data, k = 3)
mycolours_estimated <- c("darkseagreen", "orchid3", "orange")
layout(mat = cbind(c(1, 2)), heights = c(2, 1))
par(mar = c(5, 5, 1, 1))
ids <- c(589, 229, 106)
plot(simul$data,
  pch = 19, cex = 3,
  col = mycolours_estimated[pam$clustering],
  xlab = "First attribute",
  ylab = "Second attribute",
  cex.lab = 1.5, las = 1
)
text(simul$data, labels = 1:nrow(simul$data))

points(simul$data,
  pch = 19, cex = 3,
  col = mycolours_estimated[pam$clustering],
  xlab = "X1", ylab = "X2", cex.lab = 1.5, las = 1
)
text(simul$data, labels = 1:nrow(simul$data))
legend("topright",
  col = mycolours_estimated, pch = 15,
  legend = paste0("Estimated cluster ", 1:3),
  bty = "n", cex = 1.5
)

true_sil <- silhouette(x = pam$clustering, dist = dist(simul$data))[, 3]
plot(true_sil,
  type = "h", lwd = 5, lend = 1, bty = "n",
  las = 1, xlab = "", xaxt = "n", ylim = c(-0.4, 0.7),
  ylab = "Silhouette width", cex.lab = 1.5,
  col = mycolours_estimated[pam$clustering],
  panel.first = abline(v = 1:nrow(simul$data), lty = 2, col = mycolours_estimated[pam$clustering])
)
for (k in 1:nrow(simul$data)) {
  axis(
    side = 1, at = k, labels = k, las = 2, col.axis = mycolours_estimated[pam$clustering[k]],
    tick = FALSE, line = -0.5
  )
}
for (k in 1:length(unique(simul$theta))) {
  axis(side = 1, at = range(which(simul$theta == k)) + c(-0.25, 0.25), labels = NA, line = 2)
  axis(
    side = 1, at = mean(which(simul$theta == k)), labels = paste0("Simulated cluster ", k),
    line = 1.5, col.axis = darken(mycolours[k], amount = 0.5), tick = FALSE
  )
}
dev.off()

# Calibration by silhouette coefficient
score <- rep(NA, 20)
for (k in 2:20) {
  set.seed(1)
  pam <- pam(x = simul$data, k = k)
  myclusters <- pam$clustering
  mysilhouette <- silhouette(x = myclusters, dist = dist(simul$data))
  score[k] <- mean(mysilhouette[, 3])
}
plot(score)
