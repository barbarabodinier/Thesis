rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter3/")

library(sharp)
library(rCOSA)

theta_xc <- rbind(
  c(rep(1, 5), rep(0, 5)),
  c(rep(1, 4), rep(0, 6)),
  c(rep(1, 4), rep(0, 6))
)
set.seed(1)
simul <- SimulateClustering(
  n = c(20, 20, 20),
  theta_xc = theta_xc,
  ev_xc = 0.9
)
simul$theta_xc
rslt_dflt_cosa2 <- cosa2(X = data.frame(simul$data))
rslt_dflt_cosa2$D
rslt_dflt_cosa2$W

Heatmap(as.matrix(dist(simul$data)))
Heatmap(as.matrix(rslt_dflt_cosa2$D))
Heatmap(as.matrix(rslt_dflt_cosa2$W))

plot(
  as.matrix(dist(x))[upper.tri(as.matrix(dist(x)))],
  mydist[upper.tri(mydist)]
)

W <- as.matrix(rslt_dflt_cosa2$W)
x <- simul$data
mydist <- matrix(0, nrow = nrow(x), ncol = nrow(x))
par(mfrow = c(2, 5))
for (k in 1:ncol(W)) {
  maxW <- matrix(NA, nrow = nrow(W), ncol = nrow(W))
  for (i in 1:nrow(W)) {
    for (j in 1:nrow(W)) {
      maxW[i, j] <- max(W[i, k], W[j, k])
    }
  }
  tmpdist <- as.matrix(dist(x[, k, drop = FALSE]))
  mydist <- mydist + maxW * tmpdist
  Heatmap(maxW * tmpdist)
}
mydist[1:5, 1:5]

apply(W, 1, sum)
for (k in 1:ncol(W)) {
  boxplot(split(W[, k], f = simul$theta), ylim = range(W))
}

plot(
  as.matrix(rslt_dflt_cosa2$D)[upper.tri(as.matrix(rslt_dflt_cosa2$D))],
  mydist[upper.tri(mydist)]
)
