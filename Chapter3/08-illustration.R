rm(list = ls())

library(sharp)
library(igraph)
library(randomcoloR)
library(colorspace)
library(aricode)

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
stab <- Clustering(xdata = x, implementation = HierarchicalClustering)

# Initialising figure
{
  pdf("Working_figures/Illustration.pdf", width = 23, height = 20)
  layout(mat = matrix(c(
    1, 5, 9, 13,
    2, 6, 10, 14,
    3, 7, 11, 15,
    4, 8, 12, 16
  ), ncol = 4, byrow = TRUE))

  # Distance metric on different subsamples
  par(mar = c(5, 5, 1, 7))
  for (k in 1:3) {
    set.seed(k)
    s <- Resample(data = x)
    xsub <- x[s, ]
    mydist <- as.matrix(dist(xsub, upper = TRUE))
    mydist_full <- matrix(NA, ncol = nrow(x), nrow = nrow(x))
    rownames(mydist_full) <- colnames(mydist_full) <- rownames(x)
    for (i in 1:nrow(mydist_full)) {
      row_id <- rownames(mydist_full)[i]
      if (row_id %in% rownames(mydist)) {
        for (j in 1:ncol(mydist_full)) {
          col_id <- colnames(mydist_full)[j]
          if (col_id %in% colnames(mydist)) {
            mydist_full[row_id, col_id] <- mydist[row_id, col_id]
          }
        }
      }
    }
    Heatmap(mydist_full,
      legend_range = c(0, 10),
      legend = ifelse(k == 3, yes = TRUE, no = FALSE)
    )

    for (l in c(2, 3, 4)) {
      myhclust <- hclust(as.dist(mydist), method = "complete")
      members <- CoMembership(cutree(myhclust, k = l))
      members_full <- matrix(NA, ncol = nrow(x), nrow = nrow(x))
      rownames(members_full) <- colnames(members_full) <- rownames(x)
      for (i in 1:nrow(members_full)) {
        row_id <- rownames(members_full)[i]
        if (row_id %in% rownames(members)) {
          for (j in 1:ncol(members_full)) {
            col_id <- colnames(members_full)[j]
            if (col_id %in% colnames(members)) {
              members_full[row_id, col_id] <- members[row_id, col_id]
            }
          }
        }
      }
      Heatmap(members_full, col = c("royalblue", "darkred"), legend = FALSE)
    }
  }

  # Consensus matrices
  plot.new()
  for (l in c(2, 3, 4)) {
    Heatmap(ConsensusMatrix(stab, argmax_id = l))
  }
  dev.off()
}


# Calibration plot
max_N <- 10
{
  pdf("Working_figures/Illustration_calibration.pdf", width = 15, height = 5)
  par(mar = c(1, 5, 5, 1))
  plot(stab$nc[1:max_N], stab$Sc[1:max_N],
    pch = 19, col = "navy",
    panel.first = c(
      abline(h = stab$Sc, lty = 3, col = "grey"),
      abline(v = stab$nc, lty = 3, col = "grey")
    ),
    xlab = "", ylab = "Consensus Score", xaxt = "n",
    cex.lab = 1.5, las = 3
  )
  axis(side = 3, at = stab$nc[1:max_N], las = 2)
  mtext(text = "Number of clusters", side = 3, line = 3, cex = 1.5)
  dev.off()
}
