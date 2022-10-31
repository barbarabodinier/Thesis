rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter2/")

library(igraph)
library(sharp)

# Data simulation
set.seed(1)
simul <- SimulateGraphical(
  nu_within = 0.2,
  output_matrices = TRUE
)

Lambda <- LambdaSequence(lmax = 0.4, lmin = 0.06, cardinal = 5)
graph_lasso <- PenalisedGraphical(xdata = simul$data, Lambda = Lambda)

pdf("Working_figures/Graphical_lasso.pdf",
  width = 25, height = 5
)
par(mfrow = c(1, 5), mar = c(3, 3, 4, 3))
mygraph <- Graph(graph_lasso$adjacency[, , 5],
  node_label = paste0("X", 1:ncol(graph_lasso$adjacency)),
  node_colour = "dodgerblue",
  edge_colour = "black",
  satellites = TRUE
)
mylayout <- layout_with_fr(mygraph)
for (k in 1:dim(graph_lasso$adjacency)[3]) {
  mygraph <- Graph(graph_lasso$adjacency[, , k],
    node_label = paste0("X", 1:ncol(graph_lasso$adjacency)),
    node_colour = "dodgerblue",
    edge_colour = "black",
    satellites = TRUE
  )
  V(mygraph)$size <- 20
  V(mygraph)$label.cex <- 2

  plot(mygraph, layout = mylayout)
  title(
    main = eval(parse(text = paste0(
      "expression(lambda*'=",
      formatC(Lambda[k], format = "f", digits = 2),
      "')"
    ))),
    cex.main = 4
  )
  text(-1.25, 1.3, label = LETTERS[k], cex = 7)
}
dev.off()
