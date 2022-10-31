rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter3/")

library(sharp)
library(mclust)
library(MASS)
library(colorspace)
library(mnormt)

set.seed(0)
n <- c(200, 200, 200)
simul <- SimulateClustering(
  n = n,
  pk = 2,
  theta_xc = c(1, 1),
  adjacency = matrix(c(0, 1, 1, 0), byrow = TRUE, ncol = 2),
  ev_xx = 0.9,
  v_within = 1,
  v_sign = -1,
  ev_xc = 0.7,
  output_matrices = TRUE
)

pdf(paste0("Working_figures/Multivariate_gaussian_mixture_example.pdf"),
    width = 14, height = 7
)
mycolours=c("royalblue","tomato","tan")
par(mar = c(5, 5, 1, 1), mfrow=c(1,2))
plot(simul$data,
     pch = 19,
     col = NA,
     xlab = "X1", ylab = "X2", cex.lab = 1.5, las = 1
)

for (k in 1:nrow(simul$mu_mixture)){
  x     <- seq(-5, 5, length.out=1000)
  y     <- seq(-5, 5, length.out=1000)
  mu    <- simul$mu_mixture[k,]
  sigma <- simul$sigma
  f     <- function(x, y) dmnorm(cbind(x, y), mu, sigma)
  z     <- outer(x, y, f)

  contour(x, y, z,
          las = 1, cex.main = 1.5,
          nlevels = 5,
          xlab = "", ylab = "", bty = "n",
          lwd = 2, cex.lab = 1.5,
          col = darken(mycolours[k], amount=0.5),
          pch = 19, cex = 0.5, las = 1, lty = 1,
          drawlabels = TRUE, add=TRUE)
}

plot(simul$data,
     pch = 19,
     col = mycolours[simul$theta], cex=0.5,
     xlab = "X1", ylab = "X2", cex.lab = 1.5, las = 1
)
legend("bottomright",
       pch = 19, col = mycolours, bty = "n",
       cex = 1.5, pt.cex = 0.5,
       legend = paste("Cluster", 1:length(n))
)
dev.off()

# Estimation of the Gaussian mixture model
gmm=Mclust(simul$data, G=3)

pdf(paste0("Working_figures/Multivariate_gaussian_mixture_example_Mclust.pdf"),
    width = 7, height = 7
)
mycolours=c("royalblue","tomato","tan")
par(mar = c(5, 5, 1, 1))
ids=c(589, 229, 106)
plot(simul$data,
     pch = 19, cex=0.5,
     col = mycolours[simul$theta],
     xlab = "X1", ylab = "X2", cex.lab = 1.5, las = 1
)

for (k in 1:ncol(gmm$parameters$mean)){
  x     <- seq(-5, 5, length.out=1000)
  y     <- seq(-5, 5, length.out=1000)
  mu    <- gmm$parameters$mean[,k]
  sigma <- gmm$parameters$variance$sigma[,,k]
  f     <- function(x, y) dmnorm(cbind(x, y), mu, sigma)
  z     <- outer(x, y, f)

  contour(x, y, z,
          las = 1, cex.main = 1.5,
          nlevels = 5,
          xlab = "", ylab = "", bty = "n",
          lwd = 2, cex.lab = 1.5,
          col = darken(mycolours[k], amount=0.5),
          pch = 19, cex = 0.5, las = 1, lty = 2,
          drawlabels = TRUE, add=TRUE)
}

points(simul$data[ids,,drop=FALSE], cex=2, lwd=2)
text(simul$data[ids,,drop=FALSE],
     labels=c(1,2,3),
     pos=c(3,4,1), cex=1.5, offset = 0.7)
legend("bottomright",
       pch = 19, col = mycolours, bty = "n",
       cex = 1.5, pt.cex = 0.5,
       legend = paste("Cluster", 1:length(n))
)
dev.off()

gmm$z[ids,]

