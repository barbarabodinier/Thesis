rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter2/")

library(igraph)
library(sharp)

# Random network
set.seed(1)
p <- 20
simul <- SimulateAdjacency(pk = p, topology = "random", nu_within = 0.3)
mygraph1 <- Graph(simul,
  node_label = as.character(1:p),
  node_colour = "dodgerblue",
  edge_colour = "black"
)
V(mygraph1)$size <- 15
V(mygraph1)$label <- rep("", length(V(mygraph1)$label))

# Scale-free network
set.seed(1)
p <- 20
simul <- SimulateAdjacency(pk = p, topology = "scale-free")
mygraph2 <- Graph(simul,
  node_label = as.character(1:p),
  node_colour = "dodgerblue",
  edge_colour = "black"
)
V(mygraph2)$size <- 15
V(mygraph2)$label <- rep("", length(V(mygraph2)$label))

# Network with communities
set.seed(1)
p <- 20
simul <- SimulateAdjacency(
  pk = c(5, 5, 5, 5), topology = "random",
  nu_within = 0.95, nu_between = 0.05
)
mygraph3 <- Graph(simul,
  node_label = as.character(1:p),
  node_colour = c(
    rep("tomato", 5),
    rep("darkseagreen", 5),
    rep("tan", 5),
    rep("steelblue", 5)
  ),
  edge_colour = "black"
)
V(mygraph3)$size <- 15
V(mygraph3)$label <- rep("", length(V(mygraph3)$label))

# Network with connected components
set.seed(1)
p <- 20
simul <- SimulateAdjacency(
  pk = c(5, 5, 5, 5), topology = "random",
  nu_within = 0.9, nu_between = 0
)
mygraph4 <- Graph(simul,
  node_label = as.character(1:p),
  node_colour = c(
    rep("tomato", 5),
    rep("darkseagreen", 5),
    rep("tan", 5),
    rep("steelblue", 5)
  ),
  edge_colour = "black"
)
V(mygraph4)$size <- 15
V(mygraph4)$label <- rep("", length(V(mygraph4)$label))

pdf("Working_figures/Network_topology_panel.pdf",
  width = 14, height = 14
)
par(mfrow = c(2, 2), mar = rep(5, 4))
plot(mygraph1, layout = layout_with_kk(mygraph1))
plot(mygraph2, layout = layout_with_kk(mygraph2))
set.seed(1)
plot(mygraph3, layout = layout_with_fr(mygraph3))
set.seed(1)
plot(mygraph4, layout = layout_with_fr(mygraph4))
dev.off()
