rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter2/")

library(igraph)
library(sharp)

# Undirected graph
set.seed(1)
simul <- SimulateAdjacency(pk = 5, topology = "scale-free")
mygraph <- Graph(simul,
  node_label = as.character(1:5),
  node_colour = "dodgerblue",
  edge_colour = "black"
)
V(mygraph)$size <- 20
V(mygraph)$label.cex <- 2

pdf("Working_figures/Undirected_graph.pdf")
par(mar = rep(0, 4))
plot(mygraph, layout = layout_with_kk(mygraph))
dev.off()

get.adjacency(mygraph)

# Directed graph
set.seed(1)
simul <- SimulateAdjacency(pk = 5, topology = "scale-free")
simul[2, 3] <- 0
simul[2, 5] <- 0
simul[4, 2] <- 0
simul[1, 4] <- 0
mygraph <- Graph(simul,
  mode = "directed",
  node_label = as.character(1:5),
  node_colour = "dodgerblue",
  edge_colour = "black"
)
V(mygraph)$size <- 20
V(mygraph)$label.cex <- 2
E(mygraph)$arrow.size <- 1

pdf("Working_figures/Directed_graph.pdf")
par(mar = rep(0, 4))
plot(mygraph, layout = layout_with_kk(mygraph))
dev.off()

get.adjacency(mygraph)
