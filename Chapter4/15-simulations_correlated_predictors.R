rm(list = ls())

library(fake)
library(sharp)

# Reading arguments
seed <- 1
PFER_thr <- 5
pi_list <- seq(0.6, 0.9, by = 0.05)
pk <- rep(100, 10)
nu_xy <- 0.02

nu_within_list <- c(0, 0.02)
for (k in 1:length(nu_within_list)) {
  nu_within <- nu_within_list[k]
  print(nu_within)

  # Simulation
  set.seed(seed)
  xsimul <- SimulateGraphical(
    n = 500,
    pk = pk,
    nu_within = nu_within,
    nu_between = 0
  )
  simul <- SimulateRegression(
    xdata = xsimul$data,
    q = 1,
    nu_xy = nu_xy,
    beta_abs = 1,
    beta_sign = 1,
    continuous = FALSE,
    ev_xy = 0.4
  )
  # summary(lm(simul$ydata~simul$xdata[, which(simul$theta==1)]))

  # stab=VariableSelection(xdata=simul$xdata, ydata=simul$ydata)
  # SelectionPerformance(stab, simul)

  # # Simulation
  # set.seed(seed)
  # adjacency <- SimulateAdjacency(
  #   pk = pk,
  #   nu_within = nu_within,
  #   nu_between = 0
  # )
  # theta <- fake:::SamplePredictors(pk = pk, q = length(pk), nu = nu_xz, orthogonal = TRUE)
  # theta <- apply(theta, 1, sum)
  # simul <- SimulateRegression(
  #   n = 500, pk = sum(pk),
  #   adjacency_x = adjacency,
  #   theta_xz = cbind(theta),
  #   v_within = 1,
  #   eta_set = 1,
  #   family = "gaussian", ev_xz = 0.4
  # )
  print(table(simul$theta))

  # Heatmap of correlations between predictors
  mycor <- cor(simul$xdata)
  plotname <- paste0("Working_figures/Heatmap_correlation_predictors_in_regression_", k, ".pdf")
  {
    pdf(plotname, width = 8, height = 8)
    par(mar = c(5, 5, 5, 7))
    Heatmap(mycor,
      col = c("navy", "white", "darkred"),
      legend_range = c(-1, 1), legend = TRUE, axes = FALSE
    )
    dev.off()
  }
}
