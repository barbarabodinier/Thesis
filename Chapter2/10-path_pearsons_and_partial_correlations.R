rm(list = ls())

library(corpcor)
library(igraph)
library(fake)
library(sharp)

# Parameters
pk <- c(2, 3, 2)

# Preparing the directed adjacency matrix
theta <- matrix(0, nrow = sum(pk), ncol = sum(pk))
theta[1, 4] <- 1
theta[2, 6] <- 1
theta[3, 6] <- 1
theta[4, 6] <- 1
theta[2, 5] <- 1

# Data simulation
set.seed(1)
simul <- SimulateSCM(pk = pk, theta = theta, output_matrices = TRUE)

# Extracting the DAG
mygraph <- plot(simul)
mylayout <- layout_with_sugiyama(mygraph, layers = rep.int(1:length(pk), times = pk), weights = NA)
mylayout$layout[1, 1] <- 0
mylayout$layout[2, 1] <- 3
mylayout$layout[c(4, 6), 1] <- 1.5
mylayout$layout[7, 1] <- 3
mylayout$layout[5, 1] <- 3
plot.igraph(mygraph, layout = mylayout)
g_path <- mygraph
V(g_path)$color <- "dodgerblue"
V(g_path)$frame.color <- "dodgerblue"

# Matrix calculations
path_coef <- simul$path_coef
sigma <- simul$sigma
phi <- cor2pcor(simul$sigma)
V(g_path)$name <- rownames(path_coef) <- colnames(path_coef) <- rownames(sigma) <- colnames(sigma) <- rownames(phi) <- colnames(phi) <- paste0("X", 1:ncol(sigma))

# Saving figure
pdf("Working_figures/Path_coefficients_and_correlations.pdf",
  width = 18, height = 11
)
par(mfrow = c(2, 3), mar = c(5, 5, 5, 9))
Heatmap(path_coef,
  col = c("navy", "white", "red"),
  legend_range = c(-1, 1),
  text = TRUE,
  format = "f", digits = 2, zero.print = "0.00",
  cex = 1.5, cex.axis = 2,
  legend = FALSE
)
Heatmap(sigma,
  col = c("navy", "white", "red"),
  legend_range = c(-1, 1),
  text = TRUE,
  format = "f", digits = 2, zero.print = "0.00",
  cex = 1.5, cex.axis = 2,
  legend = FALSE
)
Heatmap(phi,
  col = c("navy", "white", "red"),
  legend_range = c(-1, 1),
  text = TRUE,
  format = "f", digits = 2, zero.print = "0.00",
  cex = 1.5, cex.axis = 2, cex.legend = 2
)
adjacency_sigma <- ifelse(zapsmall(sigma) != 0, yes = 1, no = 0)
diag(adjacency_sigma) <- 0
g_sigma <- Graph(adjacency_sigma,
  satellites = TRUE,
  node_colour = "dodgerblue", edge_colour = "black"
)
adjacency_phi <- ifelse(zapsmall(phi, digits = 5) != 0, yes = 1, no = 0)
diag(adjacency_phi) <- 0
g_phi <- Graph(adjacency_phi,
  satellites = TRUE,
  node_colour = "dodgerblue", edge_colour = "black"
)
V(g_path)$size <- V(g_sigma)$size <- V(g_phi)$size <- 30
V(g_path)$label.cex <- V(g_sigma)$label.cex <- V(g_phi)$label.cex <- 2
set.seed(1)
g_comp1 <- GraphComparison(g_phi, g_path, satellites = TRUE)
E(g_phi)$color[which(E(g_comp1)$color == "tomato")] <- "orange"
edgelist_phi <- get.edgelist(g_phi)
rownames(edgelist_phi) <- paste0(edgelist_phi[, 1], "-", edgelist_phi[, 2])
edgelist_sigma <- get.edgelist(g_sigma)
rownames(edgelist_sigma) <- paste0(edgelist_sigma[, 1], "-", edgelist_sigma[, 2])
E(g_sigma)$color[which(rownames(edgelist_sigma) %in% rownames(edgelist_phi)[which(E(g_comp1)$color == "tomato")])] <- "orange"
g_comp <- GraphComparison(g_sigma, g_path, satellites = TRUE)
E(g_sigma)$color[which(E(g_comp)$color == "tomato")] <- "red"
E(g_path)$color <- "black"
E(g_path)$width <- E(g_sigma)$width <- E(g_phi)$width <- 2
par(mar = c(5, 5, 2, 9))
plot(g_path, layout = mylayout)
plot(g_sigma, layout = mylayout)
plot(g_phi, layout = mylayout)
dev.off()

plotname <- "Working_figures/Equivalence_class.pdf"
pdf(plotname, width = 18, height = 11)
par(mfrow = c(2, 3), mar = c(5, 5, 5, 9))
for (k in 1:3) {
  print(k)
  # Preparing the directed adjacency matrix
  theta <- matrix(0, nrow = sum(pk), ncol = sum(pk))
  if (k == 2) {
    theta[1, 4] <- 1
  }
  if (k %in% c(1, 3)) {
    theta[4, 1] <- 1
  }
  theta[2, 6] <- 1
  theta[3, 6] <- 1
  theta[4, 6] <- 1
  if (k %in% c(1, 2)) {
    theta[2, 5] <- 1
  } else {
    theta[5, 2] <- 1
  }

  # Data simulation
  set.seed(1)
  simul <- SimulateSCM(pk = pk, theta = theta, output_matrices = TRUE)

  # Extracting the DAG
  mygraph <- Graph(
    adjacency = simul$theta,
    mode = "directed",
    node_label = paste0("X", 1:sum(pk)),
    node_colour = "dodgerblue",
    edge_colour = "black",
    satellites = TRUE
  )
  mylayout <- layout_with_sugiyama(mygraph, layers = rep.int(1:length(pk), times = pk), weights = NA)
  mylayout$layout[1, 1] <- 0
  mylayout$layout[2, 1] <- 3
  mylayout$layout[c(4, 6), 1] <- 1.5
  mylayout$layout[7, 1] <- 3
  mylayout$layout[5, 1] <- 3
  g_path <- mygraph
  V(g_path)$size <- V(g_sigma)$size <- V(g_phi)$size <- 30
  V(g_path)$label.cex <- V(g_sigma)$label.cex <- V(g_phi)$label.cex <- 2
  E(g_path)$arrow.size <- 1
  E(g_path)$width <- 2

  if (k == 1) {
    E(g_path)$color[4] <- "red"
  }
  if (k == 2) {
    E(g_path)$color[2] <- "red"
  }
  if (k == 3) {
    E(g_path)$color[3] <- "red"
    E(g_path)$color[5] <- "red"
  }

  plot(g_path, layout = mylayout)

  # adjacency_phi <- ifelse(zapsmall(phi, digits = 5) != 0, yes = 1, no = 0)
  # diag(adjacency_phi) <- 0
  # g_phi <- Graph(adjacency_phi,
  #                satellites = TRUE,
  #                node_colour = "dodgerblue", edge_colour = "black"
  # )
  # plot(g_phi, layout=mylayout)
}
dev.off()
system(paste("pdfcrop --margin 10", plotname, plotname))


plotname <- "Working_figures/Refined_dag.pdf"
pdf(plotname, width = 18, height = 11)
par(mfrow = c(2, 3), mar = c(5, 5, 5, 9))
# Preparing the directed adjacency matrix
theta <- matrix(0, nrow = sum(pk), ncol = sum(pk))
theta[3, 6] <- 1
theta[4, 6] <- 1
theta[2, 5] <- 1

# Data simulation
set.seed(1)
simul <- SimulateSCM(pk = pk, theta = theta, output_matrices = TRUE)

# Extracting the DAG
mygraph <- Graph(
  adjacency = simul$theta,
  mode = "directed",
  node_label = paste0("X", 1:sum(pk)),
  node_colour = "dodgerblue",
  edge_colour = "black",
  satellites = TRUE
)
g_path <- mygraph
V(g_path)$size <- V(g_sigma)$size <- V(g_phi)$size <- 30
V(g_path)$label.cex <- V(g_sigma)$label.cex <- V(g_phi)$label.cex <- 2
E(g_path)$arrow.size <- 1
E(g_path)$width <- 2

plot(g_path, layout = mylayout)
dev.off()
system(paste("pdfcrop --margin 10", plotname, plotname))
