rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter4/")

library(fake)
library(sharp)
library(colorspace)

# Scenario 1
for (simul_study_id in 1:3) {
  print(paste0("Simulation study ", simul_study_id))
  params_list <- read.table(paste0("Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"),
    sep = "\t", header = TRUE, stringsAsFactors = FALSE
  )

  # Heatmap of pairwise Euclidian distances
  {
    pdf(paste0("Working_figures/Heatmaps_examples_simul_study_", simul_study_id, ".pdf"),
      width = 14, height = 4
    )
    par(mfrow = c(1, nrow(params_list)), mar = c(5, 5, 2, 6))
    for (params_id in 1:nrow(params_list)) {
      # Extracting simulation parameters
      nc <- params_list[params_id, "nc"]
      equal_size <- params_list[params_id, "equal_size"]
      n_tot <- params_list[params_id, "n_tot"]
      p <- params_list[params_id, "p"]
      ev_xc <- params_list[params_id, "ev_xc"]
      nu_xc <- params_list[params_id, "nu_xc"]
      v_min <- params_list[params_id, "v_min"]
      v_max <- params_list[params_id, "v_max"]

      # Data simulation
      set.seed(1)
      if (equal_size) {
        n <- rep(1, nc) / sum(rep(1, nc)) * n_tot
      } else {
        n <- round(c(20, 50, 30, 10, 40) / sum(c(20, 50, 30, 10, 40)) * n_tot)
      }
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

      # Heatmap
      Heatmap(as.matrix(dist(simul$data)))
    }
    dev.off()
  }
}


# Scenario 2
simul_study_id <- 4
print(paste0("Simulation study ", simul_study_id))
params_list <- read.table(paste0("Simulation_parameters/Simulation_parameters_list_", simul_study_id, ".txt"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

# Heatmap of pairwise Euclidian distances
{
  pdf(paste0("Working_figures/Heatmaps_examples_simul_study_", simul_study_id, ".pdf"),
    width = 14, height = 4
  )
  par(mfrow = c(1, 3), mar = c(5, 5, 2, 6))
  params_id <- 1
  # Extracting simulation parameters
  nc <- params_list[params_id, "nc"]
  equal_size <- params_list[params_id, "equal_size"]
  n_tot <- params_list[params_id, "n_tot"]
  p <- params_list[params_id, "p"]
  ev_xc <- params_list[params_id, "ev_xc"]
  nu_xc <- params_list[params_id, "nu_xc"]
  v_min <- params_list[params_id, "v_min"]
  v_max <- params_list[params_id, "v_max"]

  # Data simulation
  set.seed(1)
  if (equal_size) {
    n <- rep(1, nc) / sum(rep(1, nc)) * n_tot
  } else {
    n <- round(c(20, 50, 30, 10, 40) / sum(c(20, 50, 30, 10, 40)) * n_tot)
  }
  pk <- round(rep(0.2, 5) * p)
  q <- round(nu_xc * p)
  theta_xc <- c(rep(1, q), rep(0, p - q))
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
    theta_xc = theta_xc,
    output_matrices = TRUE
  )
  simul$data <- scale(simul$data)

  # Heatmap (all attributes)
  Heatmap(as.matrix(dist(simul$data)))

  # Heatmap (contributing attributes)
  Heatmap(as.matrix(dist(simul$data[, which(theta_xc == 1)])))
  dev.off()
}
