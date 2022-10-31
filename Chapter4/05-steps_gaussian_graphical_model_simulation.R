rm(list = ls())

library(fake)
library(sharp)
library(igraph)
library(corpcor)

pdf("Working_figures/Gaussian_graphical_model_simulation_steps.pdf",
  width = 9, height = 12
)
par(mfrow = c(3, 2), mar = c(5, 5, 1, 6))

# Adjacency matrix
set.seed(0)
theta <- SimulateAdjacency(
  pk = 20,
  topology = "random",
  nu_within = 0.1
)

# Heatmap of adjacency matrix
rownames(theta) <- colnames(theta) <- 1:ncol(theta)
Heatmap(theta, legend = FALSE)

# Graph
mygraph <- Graph(theta,
  node_colour = "dodgerblue",
  edge_colour = "black",
  satellites = TRUE
)
V(mygraph)$size <- 20
V(mygraph)$label.cex <- 2
E(mygraph)$arrow.size <- 1
plot(mygraph)

# Precision matrix
simul <- SimulateGraphical(n = 100, theta = theta, output_matrices = TRUE)
omega <- simul$omega

# Heatmap of precision matrix
rownames(omega) <- colnames(omega) <- 1:ncol(omega)
Heatmap(omega,
  legend_range = c(-3, 3),
  col = c("navy", "white", "darkred")
)
plot.new()

# Correlation matrix
sigma <- cov2cor(solve(omega))

# Heatmap of correlation matrix
Heatmap(sigma,
  legend_range = c(-1, 1),
  col = c("navy", "white", "darkred")
)

# Heatmap of data
Heatmap(simul$data)
dev.off()
