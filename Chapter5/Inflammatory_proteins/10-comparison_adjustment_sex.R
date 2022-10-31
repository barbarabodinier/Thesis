rm(list = ls())

library(openxlsx)
library(colorspace)
library(basicPlotteR)
library(plotrix)
library(sharp)

# Loading the results
stab_full <- readRDS("Results/Inflammatory_proteins/Multivariate_lung_cancer_all.rds")
stab_women <- readRDS("Results/Inflammatory_proteins/Multivariate_lung_cancer.rds")
stab_men <- readRDS("Results/Inflammatory_proteins/Multivariate_lung_cancer_men.rds")

# Extracting selection proportions
selprop_full <- SelectionProportions(stab_full)
selprop_women <- SelectionProportions(stab_women)
selprop_men <- SelectionProportions(stab_men)
print(all(names(selprop_full) == names(selprop_women)))
print(all(names(selprop_full) == names(selprop_men)))

ids_selected_women <- unique(names(c(
  which(SelectedVariables(stab_full) == 1),
  which(SelectedVariables(stab_women) == 1)
)))
ids_selected_men <- unique(names(c(
  which(SelectedVariables(stab_full) == 1),
  which(SelectedVariables(stab_men) == 1)
)))

pdf("Working_figures/Inflammatory_proteins/Comparison_adjustment_on_sex.pdf",
  width = 14, height = 7
)
par(mfrow = c(1, 1), mar = c(5, 5, 1, 5))
plot(NA,
  xlim = c(-1, 1), ylim = c(0, 1),
  xlab = "", ylab = "", xaxt = "n",
  cex.lab = 1.5, las = 1,
  panel.first = c(
    abline(v = 0, lty = 1),
    abline(h = Argmax(stab_full)[2], lty = 2, col = "darkgreen"),
    abline(v = -Argmax(stab_men)[2], lty = 2, col = "navy"),
    abline(v = Argmax(stab_women)[2], lty = 2, col = "darkred"),
    abline(h = seq(0, 1, by = 0.2), lty = 3, col = "grey80"),
    abline(v = seq(-1, 1, by = 0.2), lty = 3, col = "grey80")
  )
)
axis(side = 1, at = seq(-1, 1, by = 0.2), labels = abs(seq(-1, 1, by = 0.2)))
axis(side = 4, at = axTicks(side = 2), labels = axTicks(side = 2), las = 1)
polygon(
  x = c(-2, -Argmax(stab_men)[2], -Argmax(stab_men)[2], -2),
  y = c(-1, -1, 2, 2), border = NA,
  col = adjustcolor(col = "royalblue", alpha.f = 0.1)
)
polygon(
  x = c(2, Argmax(stab_women)[2], Argmax(stab_women)[2], 2),
  y = c(-1, -1, 2, 2), border = NA,
  col = adjustcolor(col = "tomato", alpha.f = 0.1)
)
polygon(
  x = c(-2, -2, 2, 2),
  y = c(Argmax(stab_full)[2], 2, 2, Argmax(stab_full)[2]), border = NA,
  col = adjustcolor(col = "forestgreen", alpha.f = 0.1)
)
tmpcol <- "darkred"
ids_selected <- ids_selected_women
ids_highlight <- ids_selected
points(selprop_women, selprop_full,
  pch = 19,
  col = lighten(tmpcol, amount = 0.85)
)
points(selprop_women[ids_selected],
  selprop_full[ids_selected],
  pch = 19, col = tmpcol
)
addTextLabels(selprop_women[ids_highlight],
  selprop_full[ids_highlight],
  keepLabelsInside = TRUE,
  col.label = ifelse(ids_highlight %in% ids_selected,
    yes = darken(tmpcol, amount = 0.2),
    no = lighten(tmpcol, amount = 0.6)
  ),
  col.line = ifelse(ids_highlight %in% ids_selected,
    yes = tmpcol,
    no = lighten(tmpcol, amount = 0.85)
  ),
  col.background = NA,
  border = NA,
  cex.label = 1,
  labels = ids_highlight
)
tmpcol <- "navy"
ids_selected <- ids_selected_men
ids_highlight <- ids_selected
points(-selprop_men, selprop_full,
  pch = 19,
  col = lighten(tmpcol, amount = 0.85)
)
points(-selprop_men[ids_selected],
  selprop_full[ids_selected],
  pch = 19, col = tmpcol
)
addTextLabels(-selprop_men[ids_highlight],
  selprop_full[ids_highlight],
  keepLabelsInside = TRUE,
  col.label = ifelse(ids_highlight %in% ids_selected,
    yes = darken(tmpcol, amount = 0.2),
    no = lighten(tmpcol, amount = 0.6)
  ),
  col.line = ifelse(ids_highlight %in% ids_selected,
    yes = tmpcol,
    no = lighten(tmpcol, amount = 0.85)
  ),
  col.background = NA,
  border = NA,
  cex.label = 1,
  labels = ids_highlight
)
dev.off()
