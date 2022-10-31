rm(list = ls())

library(fake)
library(corpcor)

# Data simulation
set.seed(1)
simul1 <- SimulateGraphical(n = 100, pk = 20, nu_within = 0.1, output_matrices = TRUE)
set.seed(1)
simul2 <- SimulateGraphical(n = 500, pk = 20, nu_within = 0.1, output_matrices = TRUE)
set.seed(1)
simul3 <- SimulateGraphical(n = 1000, pk = 20, nu_within = 0.1, output_matrices = TRUE)

pdf("Working_figures/True_and_empirical_partial_correlations.pdf",
  width = 18, height = 4
)
par(mar = c(1, 1, 1, 5.5), mfrow = c(1, 4))
Heatmap(
  mat = cor2pcor(cov2cor(solve(simul1$omega))),
  col = c("navy", "white", "darkred"),
  legend_range = c(-1, 1),
  legend = FALSE
)
for (k in 1:3) {
  set.seed(1)
  simul <- eval(parse(text = paste0("simul", k)))
  Heatmap(
    mat = cor2pcor(cor(simul$data)),
    col = c("navy", "white", "darkred"),
    legend_range = c(-1, 1),
    legend = ifelse(k == 3, yes = TRUE, no = FALSE)
  )
}
dev.off()
