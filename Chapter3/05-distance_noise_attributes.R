rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter3/")

library(sharp)
library(sparcl)
library(rCOSA)
library(colorspace)

# Exporting all functions from sharp (including internal ones)
r <- unclass(lsf.str(envir = asNamespace("sharp"), all = T))
for (name in r) eval(parse(text = paste0(name, "<-sharp:::", name)))
r <- unclass(lsf.str(envir = asNamespace("fake"), all = T))
for (name in r) eval(parse(text = paste0(name, "<-fake:::", name)))

# Loading all additional functions
myfunctions <- list.files("Scripts/Functions/")
myfunctions <- myfunctions[myfunctions != "Former"]
for (k in 1:length(myfunctions)) {
  source(paste0("Scripts/Functions/", myfunctions[k]))
}

pdf(paste0("Working_figures/Distance_with_noise_attributes.pdf"),
  width = 7, height = 10
)
par(mfrow = c(3, 2), mar = c(5, 3, 1, 5))

# Simulation with contributing attributes
n <- c(20, 20, 20)
set.seed(0)
simul <- SimulateClustering(
  n = n,
  ev_xc = 0.8,
  theta_xc = matrix(c(
    1, 0, 0, 1, 1,
    0, 1, 0, 1, 1,
    0, 0, 1, 1, 1
  ), byrow = TRUE, ncol = 5)
)
x <- simul$data

# Adding noise attributes
for (p_noise in c(0, 5, 5)) {
  if (p_noise != 0) {
    simul2 <- SimulateClustering(
      n = sum(n),
      theta_xc = c(
        rep(0, p_noise)
      )
    )
    x <- cbind(x, simul2$data)
  }

  # Heatmap of simulated data
  rownames(x) <- NULL
  colnames(x) <- NULL
  Heatmap(x,
    legend_range = c(-4, 4),
    col = c("navy", "ivory", "darkred")
  )
  for (k in 1:ncol(x)) {
    axis(
      side = 1, at = k - 0.5, las = 2,
      labels = paste0("attr", k),
      col.axis = ifelse(k <= 5, yes = "darkred", no = "grey60"),
      col.ticks = ifelse(k <= 5, yes = "darkred", no = "grey60")
    )
  }
  tmpaxis <- c(0, cumsum(n))
  axis(side = 2, at = tmpaxis, labels = NA)
  axis(
    side = 2, at = apply(rbind(tmpaxis[-1], tmpaxis[-length(tmpaxis)]), 2, mean),
    labels = paste("Cluster", length(n):1), tick = FALSE
  )

  # Heatmap of distance matrix
  tmpdist <- as.matrix(dist(x))
  rownames(tmpdist) <- colnames(tmpdist) <- NULL
  Heatmap(tmpdist)
  axis(side = 1, at = tmpaxis, labels = NA)
  axis(side = 2, at = tmpaxis, labels = NA)
  axis(
    side = 1, at = apply(rbind(tmpaxis[-1], tmpaxis[-length(tmpaxis)]), 2, mean),
    labels = paste("Cluster", 1:length(n)), tick = FALSE
  )
  axis(
    side = 2, at = apply(rbind(tmpaxis[-1], tmpaxis[-length(tmpaxis)]), 2, mean),
    labels = paste("Cluster", length(n):1), tick = FALSE
  )
}
dev.off()

# Scenario with 5 noise attributes
x <- x[, 1:10]
Heatmap(as.matrix(dist(x)))

pdf(paste0("Working_figures/Scatter_plot_clusters.pdf"),
    width = 7, height = 7
)
par(mar=c(5,5,1,1))
mycolours=c("royalblue","tomato","tan")
plot(x[,c(1,2)],
     col=mycolours[simul$theta],
     pch=19,
     xlab="First attribute",
     ylab="Second attribute",
     cex.lab=1.5, las=1)
legend("bottomright",
       pch=19, col=mycolours, bty="n",
       cex=1.5, pt.cex=1,
       legend=paste("Cluster", 1:length(n)))
dev.off()

# Sparse clustering using sparcl
perm.out <- HierarchicalSparseCluster.permute(x,
  wbounds = seq(1.1, 3, by = 0.01),
  nperms = 5
)
print(perm.out)
plot(perm.out)
sparsehc <- HierarchicalSparseCluster(
  x = x,
  wbound = perm.out$bestw,
  method = "complete"
)
formatC(sparsehc$ws, format = "f", digits = 2)

pdf("Working_figures/Distance_metric_per_attribute_sparcl.pdf",
  width = 6, height = 11
)
par(mfrow = c(5, 2), mar = c(1, 8, 1, 0.25))
mydist <- matrix(0, nrow = nrow(x), ncol = nrow(x))
for (k in 1:ncol(x)) {
  tmpdist <- as.matrix(dist(x[, k, drop = FALSE]))^2
  mydist <- mydist + sparsehc$ws[k] * tmpdist
  print(max(tmpdist))
  rownames(tmpdist) <- colnames(tmpdist) <- NULL
  Heatmap(tmpdist,
    legend_range = c(0, 30),
    legend = FALSE
  )
  axis(
    side = 2, at = 30, tick = FALSE, las = 1, cex.axis = 2,
    labels = paste0(
      ifelse(k != 1, yes = "+ ", no = "= "),
      formatC(sparsehc$ws, format = "f", digits = 2)[k], " x"
    ),
    col.axis = ifelse(sparsehc$ws[k] > 0.001, yes = "darkred", no = "grey60")
  )
}
dev.off()

pdf("Working_figures/Distance_metric_per_attribute_sparcl_legend.pdf",
  width = 6, height = 11
)
par(mfrow = c(5, 2), mar = c(1, 5, 1, 4))
Heatmap(tmpdist, legend_range = c(0, 30))
dev.off()

pdf("Working_figures/Distance_metric_sparcl.pdf",
  width = 6, height = 11
)
par(mfrow = c(5, 2), mar = c(1, 5, 1, 4))
rownames(mydist) <- colnames(mydist) <- NULL
Heatmap(mydist)
dev.off()

# Weighted clustering with COSA
rslt_dflt_cosa2 <- cosa2(X = as.data.frame(x), lambda = 0.5)

pdf("Working_figures/Boxplot_weights_cosa.pdf",
  width = 10, height = 6
)
par(mar = c(5, 5, 1, 1))
colnames(rslt_dflt_cosa2$W) <- paste0("attr", 1:ncol(rslt_dflt_cosa2$W))
mycolours <- "navy"
boxplot(rslt_dflt_cosa2$W,
  range = 0, las = 2,
  ylab = "Weights", cex.lab = 1.5,
  col = lighten(mycolours, amount = 0.5),
  whiskcol = mycolours, staplecol = mycolours,
  boxcol = "white",
  medcol = mycolours,
  whisklty = 1, whisklwd = 2,
  boxwex = 0.5
)
dev.off()

pdf("Working_figures/Distance_metric_per_attribute_cosa.pdf",
  width = 6, height = 11
)
par(mfrow = c(5, 2), mar = c(1, 8, 1, 0.25))
W <- as.matrix(rslt_dflt_cosa2$W)
mydist <- matrix(0, nrow = nrow(x), ncol = nrow(x))
for (k in 1:ncol(x)) {
  maxW <- matrix(NA, nrow = nrow(W), ncol = nrow(W))
  for (i in 1:nrow(W)) {
    for (j in 1:nrow(W)) {
      maxW[i, j] <- max(W[i, k], W[j, k])
    }
  }
  tmpdist <- as.matrix(dist(x[, k, drop = FALSE]))
  mydist <- mydist + maxW * tmpdist
  print(max(maxW * tmpdist))
  rownames(tmpdist) <- colnames(tmpdist) <- NULL
  Heatmap(maxW * tmpdist,
    legend_range = c(0, 0.67),
    legend = FALSE
  )
  axis(
    side = 2, at = 30, tick = FALSE, las = 1, cex.axis = 2,
    labels = ifelse(k != 1, yes = "+", no = "="), line = 2.5
  )
}
dev.off()

pdf("Working_figures/Distance_metric_per_attribute_cosa_legend.pdf",
  width = 6, height = 11
)
par(mfrow = c(5, 2), mar = c(1, 5, 1, 4))
Heatmap(maxW * tmpdist, legend_range = c(0, 0.67))
dev.off()

pdf("Working_figures/Distance_metric_cosa.pdf",
  width = 6, height = 11
)
par(mfrow = c(5, 2), mar = c(1, 5, 1, 4))
rownames(mydist) <- colnames(mydist) <- NULL
Heatmap(mydist)
dev.off()
