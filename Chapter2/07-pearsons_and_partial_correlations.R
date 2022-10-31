rm(list = ls())

library(corpcor)
library(igraph)
library(sharp)

# Data simulation
set.seed(1)
simul <- SimulateGraphical(
  nu_within = 0.2,
  output_matrices = TRUE
)

# Matrix calculations
sigma <- simul$sigma
phi <- cor2pcor(simul$sigma)
rownames(sigma) <- colnames(sigma) <- rownames(phi) <- colnames(phi) <- paste0("X", 1:ncol(sigma))

# Saving figure
pdf("Working_figures/Correlation_partial_correlation.pdf",
  width = 14, height = 14
)
par(mfrow = c(2, 2), mar = c(5, 5, 7, 7))
Heatmap(sigma,
  col = c("navy", "white", "red"),
  text = TRUE,
  format = "f", digits = 2, zero.print = "0.00",
  cex = 0.8,
  legend = FALSE
)
Heatmap(phi,
  col = c("navy", "white", "red"),
  text = TRUE,
  format = "f", digits = 2, zero.print = "0.00",
  cex = 0.8
)
adjacency_sigma <- ifelse(sigma != 0, yes = 1, no = 0)
diag(adjacency_sigma) <- 0
g_sigma <- Graph(adjacency_sigma,
  satellites = TRUE,
  node_colour = "dodgerblue", edge_colour = "black"
)
adjacency_phi <- ifelse(round(phi, digits = 2) != 0, yes = 1, no = 0)
diag(adjacency_phi) <- 0
g_phi <- Graph(adjacency_phi,
  satellites = TRUE,
  node_colour = "dodgerblue", edge_colour = "black"
)
V(g_sigma)$size <- V(g_phi)$size <- 30
V(g_sigma)$label.cex <- V(g_phi)$label.cex <- 2
set.seed(1)
g_comp <- GraphComparison(g_sigma, g_phi, satellites = TRUE)
E(g_sigma)$color[which(E(g_comp)$color == "tomato")] <- "red"
E(g_sigma)$width <- 2
E(g_phi)$width <- 2
mylayout <- layout_with_fr(g_comp, weights = rep(0.1, length(E(g_sigma))))
par(mar = c(5, 5, 2, 7))
plot(g_sigma, layout = mylayout)
plot(g_phi, layout = mylayout)
dev.off()
