rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter4/")

library(fake)
library(corpcor)
library(colorspace)

# theta1 <- SimulateAdjacency(pk = 100, nu_within = 0.05)
# theta2 <- SimulateAdjacency(pk = 100, topology = "scale-free")
# theta3 <- SimulateAdjacency(pk = c(70, 20, 10), nu_within = 0.5, nu_between = 0)

set.seed(1)
theta2 <- SimulateAdjacency(pk = 30, nu_within = 0.5)
theta3 <- SimulateAdjacency(pk = 30, topology = "scale-free")
theta1 <- SimulateAdjacency(pk = 30, nu_within = 0.05)
set.seed(1)
theta4 <- SimulateAdjacency(pk = c(30, 10, 5), nu_within = 0.8, nu_between = 0)

pdf("Working_figures/Positive_definite_matrices.pdf",
  width = 8, height = 12
)
layout(mat = t(matrix(1:8, nrow = 2)))
par(mar = c(1, 1, 1, 6))
for (k in 1:4) {
  theta <- eval(parse(text = paste0("theta", k)))
  omega_e <- SimulatePrecision(
    theta = theta, v_within = 1, v_sign = -1, u_list = 1e-10,
    pd_strategy = "min_eigenvalue"
  )
  Heatmap(
    mat = cor2pcor(cov2cor(solve(omega_e$omega))),
    col = c("navy", "white", "darkred"),
    legend_range = c(-1, 1),
    # col = c("white", darken("darkred", amount = 0.2)),
    # legend_range = c(0, 1),
    legend = FALSE
  )
  omega_d <- SimulatePrecision(
    theta = theta, v_within = 1, v_sign = -1, u_list = 1e-10,
    pd_strategy = "diagonally_dominant"
  )
  Heatmap(
    mat = cor2pcor(cov2cor(solve(omega_d$omega))),
    col = c("navy", "white", "darkred"),
    legend_range = c(-1, 1),
    # col = c("white", darken("darkred", amount=0.2)),
    # legend_range = c(0, 1),
    legend = ifelse(k == 1, yes = TRUE, no = FALSE)
  )
}
dev.off()
