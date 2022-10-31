rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter2/")

library(igraph)
library(sharp)

# Undirected graph
set.seed(1)
p <- 10
simul <- SimulateAdjacency(pk = p, topology = "scale-free")
mygraph <- Graph(simul,
  node_label = as.character(1:p),
  node_colour = "dodgerblue",
  edge_colour = "black"
)
V(mygraph)$size <- 20
V(mygraph)$label.cex <- 2

pdf("Working_figures/Node_centrality_A.pdf")
par(mar = rep(0, 4))
plot(mygraph, layout = layout_with_kk(mygraph))
dev.off()

pdf("Working_figures/Node_centrality_B.pdf")
par(mar = c(5, 5, 1, 1))
plot(degree(mygraph), betweenness(mygraph),
  xlab = "Degree", ylab = "Betweenness centrality",
  xlim = c(0.9, max(degree(mygraph)) + 0.1),
  ylim = c(-1, max(betweenness(mygraph)) + 1),
  las = 1, cex.lab = 1.5, xaxt = "n",
  pch = 19, cex = 7, col = "dodgerblue",
  panel.first = c(
    abline(h = seq(0, 25, by = 5), col = "grey", lty = 3),
    abline(v = seq(1, 4, by = 1), col = "grey", lty = 3)
  )
)
axis(side = 1, at = seq(1, 4, by = 1))
text(degree(mygraph), betweenness(mygraph),
  labels = ifelse(degree(mygraph) > 1, yes = as.character(1:p), no = ""),
  cex = 2, col = "grey20"
)
text(
  x = 1.27, y = 2.5, labels = " { 1, 2, 5, 6, 7 }", pos = 4,
  cex = 2, col = "grey20"
)
segments(x0 = 1, x1 = 1.32, y0 = 0, y1 = 2.5, lwd = 2)
dev.off()
