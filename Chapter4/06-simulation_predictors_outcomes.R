rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter4/")

library(fake)
library(sharp)
library(colorspace)

# Simulating predictors and outcomes
set.seed(1)
simul <- SimulateRegression(
  n = 100, pk = rep(5, 3),
  nu_within = 0.8,
  nu_xz = 0.3,
  eta = matrix(c(
    1, 1, 1,
    0, 1, 1,
    0, 0, 1
  ),
  byrow = TRUE, ncol = 3
  )
)

# Conditional independence structure
theta <- simul$adjacency
mygraph <- Graph(
  adjacency = theta,
  node_label = c(
    paste0("Y", 1:3),
    paste0("Z", 1:3),
    paste0("X", 1:15)
  ),
  node_colour = c(
    c("dodgerblue", "tomato", "forestgreen"),
    rep("tan", 3),
    rep("grey", 15)
  )
)
E(mygraph)$color <- ifelse(apply(get.edgelist(mygraph), 1, FUN = function(x) {
  any(grepl("latent", x))
}), yes = "orange", no = "black")
E(mygraph)$width <- ifelse(apply(get.edgelist(mygraph), 1, FUN = function(x) {
  any(grepl("latent", x))
}), yes = 2, no = 0.5)
E(mygraph)$color <- ifelse(apply(get.edgelist(mygraph), 1, FUN = function(x) {
  any(grepl("outcome1", x))
}), yes = darken("dodgerblue", amount = 0.3), no = E(mygraph)$color)
E(mygraph)$color <- ifelse(apply(get.edgelist(mygraph), 1, FUN = function(x) {
  any(grepl("outcome2", x))
}), yes = darken("tomato", amount = 0.3), no = E(mygraph)$color)
E(mygraph)$color <- ifelse(apply(get.edgelist(mygraph), 1, FUN = function(x) {
  any(grepl("outcome3", x))
}), yes = darken("forestgreen", amount = 0.3), no = E(mygraph)$color)
V(mygraph)$size <- 10
V(mygraph)$label.cex <- 1

pdf("Working_figures/Simulation_predictors_outcomes.pdf",
  width = 7, height = 7
)
par(mar = rep(0, 4))
set.seed(1)
plot(mygraph, layout = layout_with_fr(mygraph, weights = rep(2e-4, length(E(mygraph)))))
dev.off()
