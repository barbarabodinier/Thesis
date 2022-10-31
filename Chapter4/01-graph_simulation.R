rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter4/")

library(fake)
library(sharp)
library(igraph)
devtools::load_all("/Users/barbara/Dropbox/R_packages/Simulations/fake")

pdf("Working_figures/Random_graphs.pdf",
  width = 18, height = 6
)
par(mfrow = c(1, 3), mar = rep(4, 4))
set.seed(3)
for (nu in c(0.2, 0.5, 0.8)) {
  theta <- SimulateAdjacency(
    pk = 20,
    nu_within = nu,
    topology = "random"
  )
  mygraph <- Graph(theta,
    node_colour = "dodgerblue",
    node_label = rep("", ncol(theta)),
    edge_colour = "black",
    satellites = TRUE
  )
  V(mygraph)$size <- 20
  V(mygraph)$label.cex <- 2
  E(mygraph)$arrow.size <- 1
  plot(mygraph)
  mtext(
    text = eval(parse(text = paste0("expression(nu*'=", nu, "')"))),
    side = 3, cex = 2
  )
  mtext(
    text = eval(parse(text = paste0("expression(N[E]*'=", length(E(mygraph)), "')"))),
    side = 1, cex = 2, line = 2
  )
}
dev.off()

pdf("Working_figures/Scale_free_graph.pdf",
  width = 16, height = 8
)
par(mar = rep(0, 4))
set.seed(0)
theta <- as.matrix(huge.generator(
  n = 10, d = 30,
  graph = "scale-free",
  vis = FALSE
)$theta)
mygraph <- Graph(theta,
  node_colour = "dodgerblue",
  node_label = 1:ncol(theta),
  edge_colour = "black",
  satellites = TRUE
)
V(mygraph)$size <- 8
V(mygraph)$label.cex <- 2
E(mygraph)$arrow.size <- 1
plot(mygraph, layout = layout_with_kk(mygraph), asp = 0.5)
dev.off()

pdf("Working_figures/Graph_with_communities.pdf",
  width = 10, height = 15
)
par(mfrow = c(3, 2), mar = rep(5, 4))
set.seed(1)
for (nu in c(0.2, 0.1, 0)) {
  theta <- SimulateAdjacency(
    pk = rep(5, 4),
    nu_within = 0.8,
    nu_between = nu,
    topology = "random"
  )
  rownames(theta) <- colnames(theta) <- 1:ncol(theta)
  Heatmap(theta, legend = FALSE)
  mygraph <- Graph(theta,
    node_colour = c(
      rep("tomato", 5),
      rep("darkseagreen", 5),
      rep("tan", 5),
      rep("steelblue", 5)
    ),
    edge_colour = "black",
    satellites = TRUE
  )
  V(mygraph)$size <- 20
  V(mygraph)$label.cex <- 2
  E(mygraph)$arrow.size <- 1
  plot(mygraph)
}
dev.off()
