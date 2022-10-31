rm(list = ls())

library(fake)
library(sharp)
library(colorspace)

# Simulation parameters
# method <- "hierarchical"
method <- "pam"
for (simul_study_id in c("1")) {
  print(paste0("Simulation study ", simul_study_id))

  # Template design
  mycolours <- lighten(c(
    "tan",
    "darkgreen", lighten("darkgreen", amount = 0.3),
    "orange",
    "maroon4", lighten("maroon4", amount = 0.3),
    "navy", lighten("navy", amount = 0.3),
    "darkred"
  ),
  amount = 0.3
  )
  dimensionality <- c("", "", "")

  # Saving table
  continuous_metrics <- c("rand", "ari", "jaccard")
  integer_metrics <- c("G", "time")
  binary_metrics <- "signif"
  full_names <- c(
    "G*",
    "Silhouette",
    "GAP statistic",
    "G*",
    "Delta", "PAC",
    "RCSI (PAC)", "RCSI (entropy)",
    "Consensus score"
  )
  names(full_names) <- c(
    "hclust_star", "hclust_silhouette", "hclust_gap",
    "consensus_star", "delta", "pac", "rcsi_pac", "rcsi_entropy", "consensus_score"
  )
  if (method == "hierarchical") {
    algo_names <- c("Hierarchical", "Consensus hierarchical")
  }
  if (method == "kmeans") {
    algo_names <- c("K-means", "Consensus K-means")
  }
  if (method == "pam") {
    algo_names <- c("PAM", "Consensus PAM")
  }
  n_simul <- length(list.files(path = paste0("Results/Consensus_clustering/Simulations_consensus_", method, "/Simulations_", simul_study_id)))

  # Saving figure
  for (metric in c("G", "ari")) {{ pdf(paste0("Working_figures/Boxplot_", method, "_", metric, "_", simul_study_id, ".pdf"),
    width = 12, height = 4.5
  )
  par(mar = c(8, 5, 3, 1), mfrow = c(1, 3))
  for (simul_id in 1:n_simul) {
    performances <- readRDS(paste0("Results/Consensus_clustering/Simulations_consensus_", method, "/Simulations_", simul_study_id, "/Performances_", simul_id, "_merged.rds"))
    mylist <- list()
    for (k in 1:nrow(performances)) {
      mylist <- c(mylist, list(as.numeric(performances[k, metric, ])))
      assign(paste0("median", k), median(as.numeric(performances[k, metric, ])))
    }
    xseq <- (1:length(full_names))
    if (metric == "ari") {
      ylim <- c(0, 1)
    } else {
      ylim <- c(0, 30)
    }
    myylab <- ifelse(simul_id == 1, yes = toupper(metric), no = "")
    boxplot(
      at = xseq, mylist, col = mycolours, boxcol = mycolours, whiskcol = mycolours, staplecol = mycolours, medcol = darken(mycolours, amount = 0.4),
      whisklty = 1, range = 0, las = 1, main = dimensionality[simul_id], cex.main = 1.5,
      ylab = myylab, cex.lab = 2, xaxt = "n", ylim = ylim, boxwex = 0.35
    )
    abline(h = axTicks(2), lty = 3, col = "grey")
    zseq <- c(0.5, 3.5, 9.5)
    abline(v = zseq, lty = 2, col = "black")
    boxplot(
      at = xseq, mylist, col = mycolours, boxcol = mycolours, whiskcol = mycolours, staplecol = mycolours, medcol = darken(mycolours, amount = 0.4),
      whisklty = 1, range = 0, las = 1, add = TRUE,
      ylab = myylab, cex.lab = 1.5, xaxt = "n", boxwex = 0.35
    )
    for (id in c(1, 4, 9)) {
      abline(h = eval(parse(text = paste0("median", id))), col = darken(mycolours[id], amount = 0.4), lty = 2)
    }
    axis(side = 1, at = xseq, labels = full_names[rownames(performances)], las = 2)
    axis(side = 3, at = zseq, labels = NA)
    axis(side = 3, at = apply(rbind(zseq[-1], zseq[-length(zseq)]), 2, mean), labels = algo_names, tick = FALSE, cex.axis = 1.5)
  }
  dev.off() }}
}
