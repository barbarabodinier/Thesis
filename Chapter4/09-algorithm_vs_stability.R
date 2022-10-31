rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter4/")

library(fake)
library(sharp)
library(igraph)

source("Scripts/additional_functions_specific_to_comparisons_stability.R")


# Simulation parameters
n <- 200
pk <- 100
topology <- "random"
nu <- 0.02
seed <- 3

# Stability selection parameters
PFER_thr <- 10
K <- 100

# Data simulation
set.seed(seed)
simul <- SimulateGraphical(
  n = n, pk = pk, topology = topology,
  v_within = 1, nu_within = nu,
  output_matrices = TRUE
)

# Definition of the grid of penalty parameters
Lambda <- LambdaGridGraphical(xdata = simul$data, max_density = 0.9, Lambda_cardinal = 100)
Lambda_sub <- LambdaGridGraphical(xdata = simul$data, max_density = 0.5) # to subset for stability selection

# Graphical LASSO (no stability)
system.time({out_it <- CalibrateInformationTheory(x = simul$data, Lambda = Lambda)})
perf_original <- NULL
for (k in 1:length(Lambda)) {
  print(k)
  A <- out_it$path[, , k]
  perf_original <- rbind(perf_original, data.frame(SelectionPerformance(theta = A, theta_star = simul$theta)))
}
perf_original <- perf_original[perf_original$precision + perf_original$recall > 0, ]
perf_bic <- data.frame(SelectionPerformance(theta = out_it$path[, , which.min(out_it$BIC)], theta_star = simul$theta))
perf_ebic <- data.frame(SelectionPerformance(theta = out_it$path[, , which.min(out_it$EBIC)], theta_star = simul$theta))

# Stability selection
system.time({out_visited <- GraphicalModel(
  xdata = simul$data,
  Lambda = Lambda_sub
)})
system.time({
  out_visited <- GraphicalModel(
    xdata = simul$data,
    Lambda = Lambda[Lambda >= min(Lambda_sub)]
  )
})
perf_visited <- NULL
for (i in 1:nrow(out_visited$S_2d)) {
  for (j in 1:ncol(out_visited$S_2d)) {
    A <- Adjacency(out_visited, argmax_id = rbind(c(i, j)))
    perf_visited <- rbind(
      perf_visited,
      data.frame(c(SelectionPerformance(theta = A, theta_star = simul$theta),
                   lambda = as.numeric(out_visited$Lambda[i, 1]),
                   pi = out_visited$params$pi_list[j],
                   stability = out_visited$S_2d[i, j]
      ))
    )
  }
}
perf_visited <- perf_visited[perf_visited$precision + perf_visited$recall > 0, ]
myperf <- data.frame(SelectionPerformance(theta = Adjacency(out_visited), theta_star = simul$theta))

# Extracting simulation results with true number of selected variables
simul_study_id=1
topology="random"
metric="F1_score"
mylist <- list()
for (simul_id in 1:3) {
  performances <- readRDS(paste0("Results/1-graphical_model/Simulations_", simul_study_id, "_", topology, "/Performances_true_q_", simul_id, "_merged.rds"))
  for (k in 1:nrow(performances)) {
    mylist <- c(mylist, list(as.numeric(performances[k, metric, ])))
    assign(paste0("median_", simul_id, "_", k), median(as.numeric(performances[k, metric, ])))
  }
}

# Algorithm vs stability
pdf("Working_figures/Increase_in_performance_with_stability.pdf",
    width=14, height=7)
par(mar = c(5,5,1,3), mfrow=c(1,2))
mycolours <- c("tomato", "navy")

# Precision-recall 
plot(NULL,
     las = 1, cex.main = 2, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Precision", ylab = "Recall",
     cex.lab = 1.5, las = 1
)
abline(h = axTicks(2), lty = 3, col = "lightgrey", lwd = 0.5)
abline(v = axTicks(1), lty = 3, col = "lightgrey", lwd = 0.5)
points(perf_visited$precision, perf_visited$recall,
       col = mycolours[1], pch = 18, cex = 0.7
)
points(perf_original$precision, perf_original$recall,
       col = mycolours[2], pch = 16, cex = 1
)
lines(perf_original$precision, perf_original$recall,
      col = mycolours[2], lty = 2
)
legend("bottomleft", 
       legend=c("Graphical LASSO", "Stability selection"),
       col=c("navy", "tomato"), 
       pch=c(16, 18), pt.cex=c(1, 0.7), lty=c(2, NA),
       bg="white", cex=1.3)

# Distribution of F1-score
xseq=c(1,2,4,5,7,8)
mycolours=lighten(c("navy", "tomato"), amount = 0.4)
boxplot(
  at = xseq, mylist, col = mycolours, boxcol = "white", whiskcol = mycolours, staplecol = mycolours,
  medcol = darken(mycolours, amount = 0.5),
  whisklty = 1, range = 0, las = 2, 
  ylab = expression(F[1]-score), cex.lab = 1.5, xaxt = "n", boxwex = 0.6,
  xlim=c(0,9), ylim = c(0, 1)
)
for (simul_id in 1:3){
  abline(h=eval(parse(text=paste0("median_", simul_id, "_1"))), 
         col=darken(mycolours[1], amount=0.5), lty=2)
  abline(h=eval(parse(text=paste0("median_", simul_id, "_2"))), 
         col=darken(mycolours[2], amount=0.5), lty=2)
}
abline(v=c(3,6), lty=3, col="grey")
axis(side=1, at=c(1.5,4.5,7.5), labels = c("Low", "Intermediate", "High"), cex.axis=1.5)
dev.off()

# Performance stability score
pdf("Working_figures/Stability_vs_f1_score.pdf",
    width = 14, height = 7)
par(mar = c(7, 5, 5, 7), mfrow=c(1,2))
CalibrationPlot(out_visited)
plot(NULL,
     las = 1,
     col = mycolours[1], pch = 18, cex = 0.7, cex.lab = 1.5, 
     xlim=range(perf_visited$stability, na.rm = TRUE), ylim = c(0, 1),
     xlab = "Stability score", ylab = expression(F[1] - score)
)
abline(h = axTicks(2), lty = 3, col = "lightgrey", lwd = 0.5)
abline(v = axTicks(1), lty = 3, col = "lightgrey", lwd = 0.5)
id <- which.max(perf_visited$stability)
abline(h = perf_visited$F1_score[id], col = "darkred", lty = 2)
abline(v = perf_visited$stability[id], col = "darkred", lty = 2)
points(perf_visited$stability, perf_visited$F1_score,
       las = 1,
       col = mycolours[1], pch = 18, cex = 0.7, cex.lab = 1.5, ylim = c(0, 1),
       xlab = "Stability score", ylab = expression(F[1] - score)
)
dev.off()
