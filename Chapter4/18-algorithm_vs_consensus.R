rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter4/")

library(fake)
library(sharp)
devtools::load_all("/Users/barbara/Dropbox/R_packages/Stability/sharp")

# Extracting simulation results with true number of selected variables
simul_study_id <- 1
metric <- "ari"
mylist <- list()
for (simul_id in 1:3) {
  performances <- readRDS(paste0("Results/Consensus_clustering/Simulations_consensus_hierarchical/Simulations_", simul_study_id, "/Performances_", simul_id, "_merged.rds"))
  for (k in c(1, 4)) {
    mylist <- c(mylist, list(as.numeric(performances[k, metric, ])))
    assign(paste0("median_", simul_id, "_", k), median(as.numeric(performances[k, metric, ])))
  }
}

# Algorithm vs stability
pdf("Working_figures/Increase_in_performance_with_consensus.pdf",
  width = 14, height = 7
)
par(mar = c(5, 5, 1, 3), mfrow = c(1, 2))
mycolours <- c("tomato", "navy")

# Distribution of F1-score
xseq <- c(1, 2, 4, 5, 7, 8)
mycolours <- lighten(c("navy", "tomato"), amount = 0.4)
boxplot(
  at = xseq, mylist, col = mycolours, boxcol = "white", whiskcol = mycolours, staplecol = mycolours,
  medcol = darken(mycolours, amount = 0.5),
  whisklty = 1, range = 0, las = 2,
  ylab = "ARI", cex.lab = 1.5, xaxt = "n", boxwex = 0.6,
  xlim = c(0, 9), ylim = c(0, 1)
)
for (simul_id in 1:3) {
  abline(
    h = eval(parse(text = paste0("median_", simul_id, "_1"))),
    col = darken(mycolours[1], amount = 0.5), lty = 2
  )
  abline(
    h = eval(parse(text = paste0("median_", simul_id, "_4"))),
    col = darken(mycolours[2], amount = 0.5), lty = 2
  )
}
abline(v = c(3, 6), lty = 3, col = "grey")
axis(side = 1, at = c(1.5, 4.5, 7.5), labels = c("E=0.6", "E=0.5", "E=0.4"), cex.axis = 1.5)

legend("bottomleft",
  legend = c(
    "Hierarchical",
    "Consensus"
  ),
  col = c("navy", "tomato"),
  pch = 15, bty = "n", cex = 1.3
)
dev.off()


# # Reading simulation parameters
# simul_study_id=1
# params_list <- read.table(paste0("Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"),
#                           sep = "\t", header = TRUE, stringsAsFactors = FALSE
# )
#
# # Extracting simulation parameters
# params_id=1
# nc <- params_list[params_id, "nc"]
# equal_size <- params_list[params_id, "equal_size"]
# n_tot <- params_list[params_id, "n_tot"]
# p <- params_list[params_id, "p"]
# ev_xc <- params_list[params_id, "ev_xc"]
# nu_xc <- params_list[params_id, "nu_xc"]
# v_min <- params_list[params_id, "v_min"]
# v_max <- params_list[params_id, "v_max"]
#
# # Data simulation
# set.seed(0)
# if (equal_size) {
#   n <- rep(1, nc) / sum(rep(1, nc)) * n_tot
# } else {
#   n <- round(c(20, 50, 30, 10, 40) / sum(c(20, 50, 30, 10, 40)) * n_tot)
# }
# pk <- round(rep(0.2, 5) * p)
# simul <- SimulateClustering(
#   n = n,
#   pk = pk,
#   ev_xc = ev_xc,
#   nu_within = 1,
#   nu_between = 0,
#   v_within = c(v_min, v_max),
#   v_between = 0,
#   v_sign = -1,
#   pd_strategy = "min_eigenvalue",
#   nu_xc = nu_xc,
#   output_matrices = TRUE
# )
# simul$data <- scale(simul$data)
#
# # Parameters
# nc_max <- 20
#
# out=HierarchicalClustering(xdata=simul$data, nc=1:nc_max)
# perf_hclust=rep(NA, nc_max)
# for (k in 1:nc_max){
#   perf_hclust[k]=ClusteringPerformance(theta=out$comembership[,,k], theta_star = simul)$ari
# }
#
# stab=Clustering(xdata=simul$data, nc=1:nc_max)
# perf_stab=rep(NA, nc_max)
# for (k in 1:nc_max){
#   perf_stab[k]=ClusteringPerformance(theta=Clusters(stab, argmax_id=k), theta_star = simul)$ari
# }
#
# plot(perf_hclust, ylim=c(0,1))
# points(perf_stab, col="red")
