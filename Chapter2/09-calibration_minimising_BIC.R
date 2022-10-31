rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter2/")

library(igraph)
library(sharp)

source("Scripts/additional_functions_specific_to_comparisons.R")

# Data simulation
set.seed(1)
simul <- SimulateGraphical(
  n = 100, pk = 20,
  topology = "scale-free",
  output_matrices = TRUE
)

# Graphical LASSO
set.seed(1)
Lambda <- LambdaSequence(lmax = 0.75, lmin = 0.00025, cardinal = 100)
graph_lasso <- PenalisedGraphical(xdata = simul$data, Lambda = Lambda)
apply(graph_lasso$adjacency, 3, sum)
ncol(simul$data) * (ncol(simul$data) - 1)

# Selection performances
perf <- NULL
for (k in 1:dim(graph_lasso$adjacency)[3]) {
  perf <- rbind(perf, SelectionPerformance(
    theta = graph_lasso$adjacency[, , k],
    theta_star = simul$theta
  ))
}

# Calibration based on Information Theory metrics
out <- CalibrateInformationTheory(x = simul$data, Lambda = Lambda)

pdf("Working_figures/Calibration_minimising_BIC.pdf",
  width = 10, height = 10
)
par(mar = c(5, 5, 1, 1), mfrow = c(2, 2))

# Calibration plot
plot(Lambda, out$BIC,
  pch = 19, col = "grey",
  xlab = expression(lambda),
  ylab = "BIC",
  cex.lab = 1.5,
  panel.first = c(
    abline(h = min(out$BIC), lty = 2, col = "darkred"),
    abline(v = Lambda[which.min(out$BIC)], lty = 2, col = "darkred")
  )
)
points(Lambda[which.min(out$BIC)], out$BIC[which.min(out$BIC)],
  pch = 18, col = "darkred", cex = 2
)

# Calibrated graph (by BIC)
set.seed(1)
mygraph <- SelectionPerformanceGraph(
  theta = graph_lasso$adjacency[, , which.min(out$BIC)],
  theta_star = simul,
  node_label = paste0("X", 1:ncol(graph_lasso$adjacency)),
  node_colour = "dodgerblue"
)
V(mygraph)$size <- 10
V(mygraph)$label.cex <- 0.7
plot(mygraph, layout = layout_with_fr(mygraph, weights = rep(0.1, length(E(mygraph)))))

# Precision-recall plot
par(xpd = FALSE)
plot(perf$precision, perf$recall,
  xlim = c(0, 1), ylim = c(0, 1),
  pch = 19, col = "grey",
  xlab = "Precision",
  ylab = "Recall",
  cex.lab = 1.5,
  panel.first = c(
    abline(v = perf[which.min(out$BIC), "precision"], lty = 2, col = "darkred"),
    abline(h = perf[which.min(out$BIC), "recall"], lty = 2, col = "darkred"),
    abline(v = perf[which.max(perf$F1_score), "precision"], lty = 2, col = "navy"),
    abline(h = perf[which.max(perf$F1_score), "recall"], lty = 2, col = "navy")
  )
)
points(perf[which.min(out$BIC), "precision"],
  perf[which.min(out$BIC), "recall"],
  pch = 18, col = "darkred", cex = 2
)
points(perf[which.max(perf$F1_score), "precision"],
  perf[which.max(perf$F1_score), "recall"],
  pch = 18, col = "navy", cex = 2
)

# Best graph (by F1-score)
set.seed(1)
mygraph <- SelectionPerformanceGraph(
  theta = graph_lasso$adjacency[, , which.max(perf$F1_score)],
  theta_star = simul,
  node_label = paste0("X", 1:ncol(graph_lasso$adjacency)),
  node_colour = "dodgerblue"
)
V(mygraph)$size <- 10
V(mygraph)$label.cex <- 0.7
plot(mygraph, layout = layout_with_fr(mygraph, weights = rep(0.1, length(E(mygraph)))))
dev.off()
