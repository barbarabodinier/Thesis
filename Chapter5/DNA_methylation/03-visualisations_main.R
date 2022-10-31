rm(list = ls())

library(openxlsx)
library(colorspace)
library(basicPlotteR)
library(plotrix)
library(sharp)

# Parameters
tmpcol <- "navy"
subtype <- "lung_cancer_all"

# Loading the results
tmpstab <- readRDS(paste0("Results/DNA_methylation/stability_", subtype, ".rds"))
tmpexpl <- readRDS(paste0("Results/DNA_methylation/Incremental_performances_", subtype, ".rds"))
tmpsmoking <- readRDS("Results/DNA_methylation/Univariate_packyears_all.rds")

# Loading the annotation
annot <- readRDS("Data/Original/hm450.rds")
tmpexpl$names[-1] <- paste0(annot[tmpexpl$names[-1], "alt.name"], " - ", tmpexpl$names[-1])

# Extracting associated with smoking
ids_smoking <- rownames(tmpsmoking)[which(tmpsmoking$pval < 0.05 / nrow(tmpsmoking))]

# Extracting selection proportions
ids_stable <- names(SelectedVariables(tmpstab))[which(SelectedVariables(tmpstab) == 1)]
selprop <- SelectionProportions(tmpstab)

# Explanatory performances
plotname <- paste0("Working_figures/DNA_methylation/Explanatory_performances_methylation_", subtype, ".pdf")
pdf(plotname, width = 20, height = 7.5)
layout(mat = rbind(c(1, 2)), widths = c(1, 2))
par(mar = c(20, 5, 5, 7))
plot.new()
tmpexpl$names[1] <- "Packyears"
IncrementalPlot(tmpexpl,
  # ylim = c(0.65, 0.85),
  ylab = "AUC", ylas = 0, sfrac = 0.002,
  col = c(
    "darkred",
    ifelse(gsub(".* - ", "", tmpexpl$names[-1]) %in% ids_stable,
      yes = tmpcol, no = lighten(tmpcol, amount = 0.85)
    )
  ),
  quantiles = c(0.25, 0.75), xcex.axis = 0.7, xgrid = TRUE
)
abline(h = median(tmpexpl$AUC[[sum(SelectedVariables(tmpstab) == 1) + 1]]), lty = 2, col = tmpcol)
for (k in 2:length(tmpexpl$names)) {
  axis(
    side = 3, at = k, las = 2, cex.axis = 0.7,
    labels = format(selprop[gsub(".* - ", "", tmpexpl$names)], format = "f", digits = 2)[k],
    col.axis = c(
      "darkred",
      ifelse(gsub(".* - ", "", tmpexpl$names[-1]) %in% ids_stable,
        yes = tmpcol, no = lighten(tmpcol, amount = 0.85)
      )
    )[k]
  )
}
dev.off()
system(paste("pdfcrop --margin 10", plotname, plotname))

# Preparing summary table
tmpsign <- apply(tmpexpl$Beta[[sum(SelectedVariables(tmpstab)) + 3]], 2, FUN = function(x) {
  prop.table(table(factor(sign(x), levels = c("-1", "1"))))
})[2, -c(1:4)]
mytable <- cbind(
  names(tmpsign),
  annot[names(tmpsign), "alt.name"],
  annot[names(tmpsign), "chr"],
  formatC(tmpsign, format = "f", digits = 2),
  formatC(selprop[names(tmpsign)], format = "f", digits = 2),
  formatC(tmpsmoking[names(tmpsign), "coef"], format = "f", digits = 2),
  formatC(tmpsmoking[names(tmpsign), "pval"], format = "e", digits = 2)
)
colnames(mytable) <- c("cpg", "gene", "chr", "prop_pos", "prop_sel", "beta", "pval")
write.table(mytable, paste0("Tables/DNA_methylation/Stably_selected_", subtype, ".txt"),
  row.names = FALSE, col.names = TRUE, quote = FALSE, eol = "££\n", sep = "&"
)
