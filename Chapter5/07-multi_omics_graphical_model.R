rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter5/")

library(corpcor)
library(colorspace)
library(igraph)
library(sharp)

# Loading the data
proteins <- readRDS("Data/Filtered_lung_cancer/Proteins.rds")
cpg <- readRDS("Data/Filtered_lung_cancer/DNA_methylation.rds")
metab <- readRDS("Data/Filtered_lung_cancer/Metabolomics.rds")
x <- readRDS("Data/Filtered_lung_cancer/OMICs_adjusted_conf.rds")
covars <- readRDS("Data/Filtered_lung_cancer/Covariates.rds")
annot <- readRDS("Data/Original/hm450.rds")

# Parameters
PFER_thr <- 150

# Defining potential mediators from single-OMICs analyses
mediators <- c("CDCP1", "cg00073090", "C586")

# Stability selection
stab <- GraphicalModel(
  xdata = x,
  pk = c(ncol(proteins), ncol(cpg), ncol(metab)),
  PFER_thr = PFER_thr
)

# Number of edges per block
for (k in 1:6) {
  adj_vect <- Adjacency(stab)[upper.tri(Adjacency(stab))]
  block_vect <- BlockMatrix(pk = c(ncol(proteins), ncol(cpg), ncol(metab)))[upper.tri(BlockMatrix(pk = c(ncol(proteins), ncol(cpg), ncol(metab))))]
  print(sum(adj_vect[which(block_vect == k)]))
}
BlockStructure(pk = c(ncol(proteins), ncol(cpg), ncol(metab)))

# Calibration plot
pdf(paste0("Working_figures/Calibration_multi_omics_graph_", PFER_thr, ".pdf"),
  width = 15, height = 11
)
par(mfrow = c(2, 3), mar = c(7, 7, 10, 7))
CalibrationPlot(stab)
dev.off()

# Defining node colours
mycolours <- lighten(c(
  rep("darkolivegreen", ncol(proteins)),
  rep("royalblue", ncol(cpg)),
  rep("tan", ncol(metab))
), amount = 0.5)
names(mycolours) <- colnames(x)
mycolours[mediators] <- "tomato"

# Defining node labels
mylabels <- ifelse(!is.na(annot[colnames(x), "alt.name"]),
  yes = paste0(annot[colnames(x), "alt.name"], " \n ", colnames(x)),
  no = colnames(x)
)

# Creating igraph object
mygraph <- Graph(stab, node_colour = mycolours, node_label = mylabels, satellites = FALSE)
V(mygraph)$size <- rep(7, length(V(mygraph)))
V(mygraph)$label.cex <- rep(0.4, length(V(mygraph)))
E(mygraph)$color <- "grey30"

pdf(paste0("Working_figures/Multi_omics_graph_", PFER_thr, ".pdf"),
  width = 8, height = 8
)
par(mfrow = c(1, 1), mar = rep(0, 4))
set.seed(1)
plot(mygraph, layout = layout_with_fr(mygraph, weights = rep(3, length(E(mygraph)))))
dev.off()

plotname <- "Working_figures/Multi_omics_graph_legend.pdf"
pdf(plotname,
  width = 8, height = 8
)
par(mfrow = c(1, 1), mar = rep(0, 4))
plot.new()
legend("top",
  pch = 19, col = c(lighten(c(
    "darkolivegreen",
    "royalblue",
    "tan"
  ), amount = 0.5), "tomato"),
  legend = c("Inflammatory protein", "CpG site", "Metabolomics", "Potential mediator"),
  bty = "n"
)
dev.off()
system(paste("pdfcrop --margin 10", plotname, plotname))

# Partial correlations
mypcor <- cor2pcor(cor(x))
Heatmap(mypcor, col = c("navy", "white", "darkred"), legend_range = c(-0.5, 0.5))
mypcor_vect <- mypcor[upper.tri(mypcor)]
mypcor_vect[which(adj_vect == 1)]
