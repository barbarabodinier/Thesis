rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter4/")

library(fake)
library(igraph)

### Simulations

# Simulation parameters
n <- 200
pk <- 100
nu <- 0.02
seed <- 3

pdf("Working_figures/Examples_graphical_simulation_scenarios.pdf",
  width = 14, height = 7
)
par(mfrow = c(1, 2), mar = rep(2, 4))

# Random network simulation
set.seed(seed)
simul <- SimulateGraphical(
  n = n, pk = pk, topology = "random",
  v_within = 1, nu_within = nu,
  output_matrices = TRUE
)

mygraph <- Graph(simul$theta,
  node_colour = rep("dodgerblue", ncol(simul$theta)),
  edge_colour = "black",
  satellites = TRUE
)
V(mygraph)$label <- rep("", length(V(mygraph)))

set.seed(1)
plot(mygraph,
  layout = layout_with_fr(mygraph,
    weights = rep(0.5, length(E(mygraph)))
  )
)

# Scale-free network simulation
set.seed(seed)
simul <- SimulateGraphical(
  n = n, pk = pk, topology = "scale-free",
  v_within = 1,
  output_matrices = TRUE
)

mygraph <- Graph(simul$theta,
  node_colour = rep("dodgerblue", ncol(simul$theta)),
  edge_colour = "black",
  satellites = TRUE
)
V(mygraph)$label <- rep("", length(V(mygraph)))

set.seed(1)
plot(mygraph,
  layout = layout_with_fr(mygraph,
    weights = rep(0.5, length(E(mygraph)))
  )
)
dev.off()
