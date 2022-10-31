rm(list = ls())

library(fake)
library(sharp)
library(igraph)
library(randomcoloR)
library(colorspace)
library(aricode)
library(FactoMineR)
library(diceR)
library(ConsensusClusterPlus)
library(M3C)
library(abind)

source("Scripts/additional_functions_specific_to_comparisons_consensus.R")

# Model parameters
K <- 100
iters <- 25 # default 25, recommended 5-100
seed <- 4

# Calibration curves
for (simul_study_id in 1) {
  print(paste0("Simulation study ", simul_study_id))
  params_list <- read.table(paste0("Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"),
    sep = "\t", header = TRUE, stringsAsFactors = FALSE
  )
  params_id <- 1
  n_tot <- params_list[params_id, "n_tot"]
  p <- params_list[params_id, "p"]
  ev_xc <- params_list[params_id, "ev_xc"]
  nu_xc <- params_list[params_id, "nu_xc"]
  v_min <- params_list[params_id, "v_min"]
  v_max <- params_list[params_id, "v_max"]

  {
    pdf(paste0("Working_figures/Calibration_curves_simul_study_", simul_study_id, ".pdf"),
      width = 14, height = 20
    )
    layout(matrix(seq(1, 5 * 3), byrow = FALSE, ncol = 3))
    par(mar = c(5, 5, 1, 6))
    for (ev_xc in seq(0.6, 0.4, by = -0.1)) {
      # Data simulation
      set.seed(seed)
      n <- round(c(20, 50, 30, 10, 40) / sum(c(20, 50, 30, 10, 40)) * n_tot)
      pk <- round(rep(0.2, 5) * p)
      simul <- SimulateClustering(
        n = n,
        pk = pk,
        ev_xc = ev_xc,
        nu_within = 1,
        nu_between = 0,
        v_within = c(v_min, v_max),
        v_between = 0,
        v_sign = -1,
        pd_strategy = "min_eigenvalue",
        nu_xc = nu_xc,
        output_matrices = TRUE
      )
      simul$data <- scale(simul$data)

      # Heatmap of pairwise distances
      Heatmap(as.matrix(dist(simul$data)))

      # Consensus clustering
      stab <- Clustering(
        xdata = simul$data,
        implementation = HierarchicalClustering,
        K = K,
        verbose = FALSE,
        nc = 1:20,
      )

      # Delta score
      delta <- DeltaAreaCDF(stab)
      id <- ManualArgmaxId(delta)
      ManualCalibPlot(delta, ylab = expression(Delta))

      # PAC score
      pac <- PAC(stab)
      id <- ManualArgmaxId(-pac)
      ManualCalibPlot(-pac, ylab = "- PAC")

      # M3C score (PAC)
      print(system.time({
        scores <- MonteCarloScore(x = simul$data, stab, objective = "PAC", iters = iters)
      }))
      rcsi_pac <- scores$RCSI # criterion to define assignments in their code
      ManualCalibPlot(rcsi_pac, ylab = "RCSI (PAC)")

      # Consensus score
      ManualCalibPlot(stab$Sc, ylab = "Consensus score")
    }
    dev.off()
  }
}


# Scatter plot with performance
for (simul_study_id in 1) {
  print(paste0("Simulation study ", simul_study_id))
  params_list <- read.table(paste0("Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"),
    sep = "\t", header = TRUE, stringsAsFactors = FALSE
  )
  params_id <- 1
  n_tot <- params_list[params_id, "n_tot"]
  p <- params_list[params_id, "p"]
  ev_xc <- params_list[params_id, "ev_xc"]
  nu_xc <- params_list[params_id, "nu_xc"]
  v_min <- params_list[params_id, "v_min"]
  v_max <- params_list[params_id, "v_max"]

  {
    pdf(paste0("Working_figures/Clustering_performance_simul_study_", simul_study_id, ".pdf"),
      width = 14, height = 20
    )
    layout(matrix(seq(1, 5 * 3), byrow = FALSE, ncol = 3))
    par(mar = c(5, 5, 1, 6))
    for (ev_xc in seq(0.6, 0.4, by = -0.1)) {
      # Data simulation
      set.seed(seed)
      n <- round(c(20, 50, 30, 10, 40) / sum(c(20, 50, 30, 10, 40)) * n_tot)
      pk <- round(rep(0.2, 5) * p)
      simul <- SimulateClustering(
        n = n,
        pk = pk,
        ev_xc = ev_xc,
        nu_within = 1,
        nu_between = 0,
        v_within = c(v_min, v_max),
        v_between = 0,
        v_sign = -1,
        pd_strategy = "min_eigenvalue",
        nu_xc = nu_xc,
        output_matrices = TRUE
      )
      simul$data <- scale(simul$data)

      # Heatmap of pairwise distances
      Heatmap(as.matrix(dist(simul$data)))

      # Consensus clustering
      stab <- Clustering(
        xdata = simul$data,
        implementation = HierarchicalClustering,
        K = K,
        verbose = FALSE,
        nc = 1:20,
      )

      # Clustering performances for different numbers of clusters
      perf <- AllPerf(stab)

      # Delta score
      delta <- DeltaAreaCDF(stab)
      id <- ManualArgmaxId(delta)
      ScatterPerf(x = delta, perf = perf, xlab = expression(Delta))

      # PAC score
      pac <- PAC(stab)
      id <- ManualArgmaxId(-pac)
      ScatterPerf(x = -pac, perf = perf, xlab = "- PAC")

      # M3C score (PAC)
      system.time({
        scores <- MonteCarloScore(x = simul$data, stab, objective = "PAC", iters = iters)
      })
      rcsi_pac <- scores$RCSI # criterion to define assignments in their code
      ScatterPerf(x = rcsi_pac, perf = perf, xlab = "RCSI (PAC)")

      # Consensus score
      ScatterPerf(x = stab$Sc, perf = perf, xlab = "Consensus score")
    }
    dev.off()
  }
}
