rm(list = ls())

library(sharp)
library(M3C)

source("Scripts/additional_functions_specific_to_comparisons.R")

# Simulation of data with clusters
set.seed(0)
n <- c(20, 50, 30)
simul <- SimulateClustering(
  n = n,
  pk = 5,
  nu_xc = 1,
  ev_xc = 0.6
)
x <- simul$data

# Consensus clustering
stab <- Clustering(xdata = simul$data)
mc <- MonteCarloScore(x = simul$data, stab)

pdf("Working_figures/Calibration_scores.pdf",
  width = 9, height = 12
)
par(mfrow = c(3, 2), mar = c(5, 5, 1, 4.5))

# Heatmap of the consensus matrix
M <- ConsensusMatrix(stab)
Heatmap(M)

# CDF and area under the CDF
thr_list <- seq(0.01, 1, by = 0.01)
CDF <- AreaUnderCDF(M, thr_list = thr_list)$CDF
plot(thr_list, CDF,
  las = 1, pch = 19, col = "navy",
  xlab = "Co-membership proportion",
  ylab = "CDF", cex.lab = 1.5, ylim = c(0, 1)
)
for (k in 2:length(CDF)) {
  polygon(
    x = c(thr_list[k - 1], thr_list[k], thr_list[k], thr_list[k - 1]),
    y = c(0, 0, CDF[k], CDF[k]),
    col = adjustcolor("navy", alpha.f = 0.5),
    border = NA
  )
}

# Area under CDF curve
area <- rep(NA, nrow(stab$nc))
for (k in stab$nc) {
  M <- ConsensusMatrix(stab, argmax_id = k)
  area[k] <- AreaUnderCDF(M)$area
}
plot(stab$nc[-1], area[-1],
  type = "h", lend = 1, lwd = 5,
  col = "navy",
  xlab = "Number of clusters",
  ylab = "Area under CDF curve",
  cex.lab = 1.5, ylim = c(0, 1),
  xaxt = "n", las = 1
)
axis(side = 1, at = stab$nc[-1])

# Delta score
delta <- DeltaArea(areas = area[-1])
plot(stab$nc[-1], delta,
  col = "navy", pch = 19,
  xlab = "Number of clusters",
  ylab = "Delta score",
  cex.lab = 1.5,
  xaxt = "n", las = 1
)
axis(side = 1, at = stab$nc[-1])

# PAC score
pac <- PAC(stab)
plot(stab$nc[-1], pac[-1],
  col = "navy", pch = 19,
  xlab = "Number of clusters",
  ylab = "PAC score",
  cex.lab = 1.5,
  xaxt = "n", las = 1
)
axis(side = 1, at = stab$nc[-1])

# Monte Carlo score
plot(stab$nc[-1], mc[-1, "RCSI"],
  col = "navy", pch = 19,
  xlab = "Number of clusters",
  ylab = "RCSI (PAC) score",
  cex.lab = 1.5,
  xaxt = "n", las = 1
)
axis(side = 1, at = stab$nc[-1])
dev.off()

pdf("Working_figures/Consensus_score.pdf",
  width = 9, height = 12
)
par(mfrow = c(3, 2), mar = c(5, 5, 1, 4.5))
# Consensus score
plot(stab$nc[-1], stab$Sc[-1],
  col = "navy", pch = 19,
  xlab = "Number of clusters",
  ylab = "Consensus score",
  cex.lab = 1.5, ylim = c(0, 1),
  xaxt = "n", las = 1
)
axis(side = 1, at = stab$nc[-1])
dev.off()
