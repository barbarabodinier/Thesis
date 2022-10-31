rm(list = ls())

library(openxlsx)
library(colorspace)
library(basicPlotteR)
library(sharp)

# Parameters
summary_method <- "centroid"
var_of_interest <- "packyears"

# Loading the results
annot <- readRDS("/Users/barbara/Dropbox/Metabolomics_lung_cancer/Data/Annotation_representative_features.rds")
stab <- readRDS("/Users/barbara/Dropbox/Metabolomics_lung_cancer/Results/2-clustering/consensus_clustering_very_fine_both.rds")
results_summarised <- readRDS(paste0("Results/Metabolomics/Univariate_", var_of_interest, ".rds"))
results_all <- readRDS(paste0("Results/Metabolomics/Univariate_", var_of_interest, "_all_features.rds"))
description <- read.xlsx(paste0("/Users/barbara/Dropbox/Metabolomics_lung_cancer/Tables/2-clustering/Description_summary_", summary_method, ".xlsx"))
rownames(description) <- description[, 1]

# Extracting the cluster membership
myclusters <- Clusters(stab)

# Extracting the minimum and maximum p-values per cluster
min_pval <- max_pval <- rep(NA, length(unique(myclusters)))
names(min_pval) <- names(max_pval) <- rownames(results_summarised)
for (i in 1:length(unique(myclusters))) {
  tmp_ids <- names(myclusters)[myclusters == i]
  min_pval[i] <- min(results_all[tmp_ids, "pval"])
  max_pval[i] <- max(results_all[tmp_ids, "pval"])
}

{
  pdf(paste0("Working_figures/Metabolomics/Scatter_plot_features_vs_clusters_", var_of_interest, "_pval.pdf"),
    width = 7, height = 7
  )
  tmpcol <- "navy"
  amount_light <- 0.85
  par(mar = c(5, 5, 1, 1))
  plot(-log10(min_pval),
    -log10(results_summarised$pval),
    pch = 20,
    col = ifelse(results_summarised$pval < 0.05 / length(unique(myclusters)),
      yes = tmpcol, no = lighten(tmpcol, amount = amount_light)
    ),
    cex = 1,
    xlab = expression(-log[10] * "(p-value) (features)"),
    ylab = expression(-log[10] * "(p-value) (clusters)"),
    cex.lab = 1.5, las = 1,
    xlim = c(0, -log10(min(min_pval, min(results_summarised$pval)))),
    ylim = c(0, -log10(min(min_pval, min(results_summarised$pval)))),
    panel.first = c(
      abline(0, 1, lty = 3),
      abline(h = -log10(0.05 / length(unique(myclusters))), lty = 2, col = "darkred"),
      abline(v = -log10(0.05 / length(unique(myclusters))), lty = 2, col = "darkred"),
      abline(v = -log10(0.05 / nrow(results_all)), lty = 2, col = "tan")
    )
  )
  points(-log10(max_pval),
    -log10(results_summarised$pval),
    pch = 20,
    col = ifelse(results_summarised$pval < 0.05 / length(unique(myclusters)),
      yes = tmpcol, no = lighten(tmpcol, amount = amount_light)
    ),
    cex = 1
  )
  for (i in 1:length(unique(myclusters))) {
    tmp_ids <- names(myclusters)[myclusters == i]
    points(-log10(results_all[tmp_ids, "pval"]),
      -log10(rep(results_summarised$pval[i], length(tmp_ids))),
      pch = 20,
      col = ifelse(results_summarised$pval[i] < 0.05 / length(unique(myclusters)),
        yes = tmpcol, no = lighten(tmpcol, amount = amount_light)
      ),
      cex = 0.5
    )
    segments(
      x0 = -log10(max_pval[i]),
      x1 = -log10(min_pval[i]),
      y0 = -log10(results_summarised$pval[i]),
      y1 = -log10(results_summarised$pval[i]),
      col = ifelse(results_summarised$pval[i] < 0.05 / length(unique(myclusters)),
        yes = tmpcol, no = lighten(tmpcol, amount = amount_light)
      )
    )
  }
  dev.off()
}

# Extracting summary statistics
ids <- which(results_summarised$pval < 0.05 / length(unique(myclusters)))
time_range <- mass_range <- rep(NA, length(ids))
for (i in 1:length(ids)) {
  tmpnames <- names(myclusters)[myclusters == ids[i]]
  tmptime <- diff(range(as.numeric(gsub(".*@", "", tmpnames))))
  tmpmass <- diff(range(as.numeric(gsub("@.*", "", tmpnames))))
  time_range[i] <- formatC(tmptime, format = "f", digits = 2)
  mass_range[i] <- formatC(tmpmass, format = "f", digits = 2)
}
out <- cbind(
  paste0("C", ids),
  table(myclusters)[ids],
  annot[paste0("C", ids)],
  time_range,
  mass_range,
  formatC(results_summarised[ids, "coef"], format = "f", digits = 2),
  formatC(results_summarised[ids, "pval"], format = "e", digits = 2),
  formatC(min_pval[ids], format = "e", digits = 2),
  formatC(max_pval[ids], format = "e", digits = 2)
)
out <- out[sort.list(as.numeric(out[, 7]), decreasing = FALSE), ]
colnames(out) <- c("ID", "N", "Representative", "Time", "Mass", "beta", "p", "Min", "Max")

write.table(out, paste0("Tables/Metabolomics/Univariate_packyears_significant.txt"),
  row.names = FALSE, col.names = TRUE, quote = FALSE, eol = "££\n", sep = "&"
)

# Focusing on clusters 88 and 243
members <- cbind(
  do.call(rbind, strsplit(names(myclusters)[which(myclusters == "88")], split = "@")),
  rbind(
    do.call(rbind, strsplit(names(myclusters)[which(myclusters == "243")], split = "@")),
    matrix("", nrow = sum(myclusters == "88") - sum(myclusters == "243"), ncol = 2)
  )
)
members[, 1] <- format(as.numeric(members[, 1]), format = "f", digits = 4)
members[, 2] <- format(as.numeric(members[, 2]), format = "f", digits = 4)
members[, 3] <- format(as.numeric(members[, 3]), format = "f", digits = 4)
members[, 4] <- format(as.numeric(members[, 4]), format = "f", digits = 4)
members[which(is.na(as.numeric(members)))] <- ""
write.table(members, paste0("Tables/Metabolomics/Members_clusters_88_243.txt"),
  row.names = FALSE, col.names = TRUE, quote = FALSE, eol = "££\n", sep = "&"
)

# Extracting the minimum and maximum beta coefficients per cluster
min_pval <- max_pval <- rep(NA, length(unique(myclusters)))
names(min_pval) <- names(max_pval) <- rownames(results_summarised)
for (i in 1:length(unique(myclusters))) {
  tmp_ids <- names(myclusters)[myclusters == i]
  min_pval[i] <- min(results_all[tmp_ids, "coef"])
  max_pval[i] <- max(results_all[tmp_ids, "coef"])
}

pdf(paste0("Working_figures/Metabolomics/Scatter_plot_features_vs_clusters_", var_of_interest, "_beta.pdf"),
  width = 7, height = 7
)
tmpcol <- "navy"
amount_light <- 0.85
par(mar = c(5, 5, 1, 1))
plot(min_pval,
  results_summarised$coef,
  pch = 20,
  col = ifelse(results_summarised$pval < 0.05 / length(unique(myclusters)),
    yes = tmpcol, no = lighten(tmpcol, amount = amount_light)
  ),
  cex = 1,
  xlab = expression(beta * " (features)"),
  ylab = expression(beta * " (clusters)"),
  cex.lab = 1.5, las = 1,
  xlim = abs(max(max_pval, max(results_summarised$coef))) * c(-1, 1),
  ylim = abs(max(max_pval, max(results_summarised$coef))) * c(-1, 1),
  panel.first = c(
    abline(h = 0, lty = 2),
    abline(v = 0, lty = 2),
    abline(0, 1, lty = 3),
    abline(h = -log10(0.05 / length(unique(myclusters))), lty = 2, col = "darkred"),
    abline(v = -log10(0.05 / length(unique(myclusters))), lty = 2, col = "darkred")
  )
)
points(max_pval,
  results_summarised$coef,
  pch = 20,
  col = ifelse(results_summarised$pval < 0.05 / length(unique(myclusters)),
    yes = tmpcol, no = lighten(tmpcol, amount = amount_light)
  ),
  cex = 1
)
for (i in 1:length(unique(myclusters))) {
  tmp_ids <- names(myclusters)[myclusters == i]
  points(results_all[tmp_ids, "coef"],
    rep(results_summarised$coef[i], length(tmp_ids)),
    pch = 20,
    col = ifelse(results_summarised$pval[i] < 0.05 / length(unique(myclusters)),
      yes = tmpcol, no = lighten(tmpcol, amount = amount_light)
    ),
    cex = 0.5
  )
  segments(
    x0 = max_pval[i],
    x1 = min_pval[i],
    y0 = results_summarised$coef[i],
    y1 = results_summarised$coef[i],
    col = ifelse(results_summarised$pval[i] < 0.05 / length(unique(myclusters)),
      yes = tmpcol, no = lighten(tmpcol, amount = amount_light)
    )
  )
}
dev.off()
