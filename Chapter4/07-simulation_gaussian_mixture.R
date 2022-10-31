rm(list = ls())

library(fake)
library(sharp)
library(colorspace)

# Making figure
mycolours <- darken(c("navy", "red", "gold"), amount = 0.2)
{
  pdf("Working_figures/Visualisation_simulations_ev.pdf",
    width = 12, height = 7
  )
  layout(mat = matrix(1:6, byrow = FALSE, nrow = 2))
  par(mar = c(5, 5, 1, 5))
  ev_list <- c(0.5, 0.7, 0.9)
  for (k in 1:3) {
    # Data simulation
    set.seed(1)
    n <- c(20, 50, 30)
    simul <- SimulateClustering(
      n = n,
      pk = 2,
      nu_xc = 1,
      ev_xc = ev_list[k]
    )

    # Heatmap of pairwise distances
    Heatmap(as.matrix(dist(scale(simul$data))),
      legend_range = c(0, 5),
      legend = ifelse(k == 3, yes = TRUE, no = FALSE)
    )

    # Scatter plot
    xcomp <- 1
    ycomp <- 2
    plot(NA,
      xlim = range(simul$data[, 1]),
      ylim = range(simul$data[, 2]),
      las = 1, cex.lab = 1.5, cex.lab = 1.5,
      xlab = paste0("Variable ", xcomp),
      ylab = paste0("Variable ", ycomp),
    )
    abline(v = axTicks(1), lty = 3, col = "grey")
    abline(h = axTicks(2), lty = 3, col = "grey")
    text(simul$data[, xcomp], simul$data[, ycomp],
      pch = 19, cex = 0.7, las = 1,
      labels = gsub("obs", "", rownames(simul$data)),
      col = mycolours[simul$theta]
    )
  }
  dev.off()
}
