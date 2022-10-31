rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter2/")

library(sharp)
library(igraph)
library(colorspace)
library(RColorBrewer)
library(data.table)


### Loading the data

# OMICS
cpg <- readRDS("~/Dropbox/Stability_selection/Data/NOWAC_MTT_smoking_159.rds")
ttx <- readRDS("~/Dropbox/Stability_selection/Data/NOWAC_TTX_smoking_208.rds")
ids <- intersect(rownames(cpg), rownames(ttx))
cpg <- cpg[ids, ]
ttx <- ttx[ids, ]
omic <- cbind(cpg, ttx)
pk <- c(ncol(cpg), ncol(ttx))

# Annotation
ttxannot <- data.frame(readRDS("~/Dropbox/Stability_selection/Data/NOWAC_TTX_annot_208.rds"))
cpgannot <- readRDS("~/Dropbox/Stability_selection/Data/NOWAC_MTT_annot_159.rds")
cpgannot <- cbind(cpgannot[, c("alt.name", "chr"), drop = FALSE], name = rownames(cpgannot))
ttxannot <- ttxannot[, c(3, 4, 2)]
colnames(ttxannot) <- c("alt.name", "chr", "name")
omicannot <- rbind(cpgannot, ttxannot)

# Heatmap of correlations
mycor <- cor(omic)
plotname <- "Working_figures/Multi_block_heatmap.pdf"
{
  pdf(plotname, width = 8, height = 8)
  par(mar = c(5, 5, 5, 7))
  Heatmap(mycor,
    col = c("darkblue", "white", "firebrick3"),
    legend_range = c(-1, 1), legend_length = 100, legend = TRUE, axes = FALSE
  )
  axis(side = 1, at = c(0, ncol(cpg)), labels = NA)
  axis(side = 1, at = mean(c(0, ncol(cpg))), labels = "DNA methylation", tick = FALSE, cex.axis = 1.5)
  axis(side = 1, at = c(ncol(cpg), ncol(omic)), labels = NA)
  axis(side = 1, at = mean(c(ncol(cpg), ncol(omic))), labels = "Gene expression", tick = FALSE, cex.axis = 1.5)
  axis(side = 2, at = c(0, ncol(ttx)), labels = NA)
  axis(side = 2, at = mean(c(0, ncol(ttx))), labels = "Gene expression", tick = FALSE, cex.axis = 1.5)
  axis(side = 2, at = c(ncol(ttx), ncol(omic)), labels = NA)
  axis(side = 2, at = mean(c(ncol(ttx), ncol(omic))), labels = "DNA methylation", tick = FALSE, cex.axis = 1.5)
  dev.off()
}
system(paste("pdfcrop --margin 10", plotname, plotname))

# Single-block stability selection
stab <- GraphicalModel(xdata = omic, max_density = 0.2, Lambda_cardinal = 30, PFER_thr = 150)
CalibrationPlot(stab)
adjacency <- Adjacency(stab)

# Define colours
colors <- lighten(colorRampPalette(brewer.pal(12, name = "Paired"))(23), amount = 0.4)
chr_number <- as.character(omicannot$chr)
chr_number[chr_number == "X"] <- 23
chr_number[chr_number == "Y"] <- 24
chr_number <- as.numeric(chr_number)

# Make igraph object
node_label <- paste(omicannot$alt.name, "\n", omicannot$name)
mygraph <- Graph(
  adjacency = adjacency, node_colour = c(rep("skyblue", ncol(cpg)), rep("lightsalmon", ncol(ttx))),
  node_label = node_label, node_shape = c(rep("square", ncol(cpg)), rep("circle", ncol(ttx)))
)
V(mygraph)$size <- 0.7 * V(mygraph)$size

# Saving figure
myasp <- 0.7
myseed <- 1
{
  pdf(paste0("Working_figures/Single_block.pdf"),
    width = 14, height = myasp * 14
  )
  par(mar = rep(0, 4))
  set.seed(myseed)
  plot(mygraph, layout = layout_with_fr(mygraph), asp = myasp)
  dev.off()
}
