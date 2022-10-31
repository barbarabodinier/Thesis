rm(list = ls())

library(glmnet)
library(sharp)

set.seed(1)
simul <- SimulateRegression(n = 200, pk = 100, nu_xz = 0.05)

lasso <- cv.glmnet(x = simul$xdata, y = simul$ydata, family = "gaussian")

pdf("Working_figures/CV_lasso.pdf")
par(mar = c(5, 5, 1, 1))
plot(lasso, cex.lab = 1.5, las = 1, col = "navy")
dev.off()
