rm(list = ls())

library(sharp)
library(igraph)
library(randomcoloR)
library(colorspace)
library(aricode)

# Simulation of data with clusters
set.seed(0)
n <- c(20, 50, 30)
simul <- SimulateClustering(
  n = n,
  pk = 5,
  nu_xc = 1,
  ev_xc = 0.6
)
x <- simul$data

# Hierarchical clustering
mydist <- as.matrix(dist(x))
myhclust <- hclust(dist(x), method = "complete")
hclustorder <- myhclust$order

{
  pdf("Working_figures/Hierarchical_clustering_dendrogram.pdf",
    width = 15, height = 6
  )
  par(mar = c(1, 1, 1, 5))
  plot(as.dendrogram(myhclust),
    xaxt = "n", leaflab = "none", axes = FALSE,
    xlab = "", sub = "", main = "", ylab = ""
  )
  axis(side = 4, at = seq(0, max(myhclust$height), by = 2), labels = NA)
  par(xpd = TRUE)
  text(length(myhclust$order) * 1.06,
    seq(0, max(myhclust$height), by = 2),
    seq(0, max(myhclust$height), by = 2),
    srt = 270
  )
  text(length(myhclust$order) * 1.09,
    mean(c(0, max(myhclust$height))),
    adj = 0.5,
    "Euclidian distance", srt = 270, cex = 1.5
  )
  rect.hclust(myhclust, k = 3, border = "navy")
  dev.off()
}

# Consensus clustering
stab <- Clustering(xdata = x, implementation = HierarchicalClustering)
plot(stab$Sc)
hat_N <- which.max(stab$Sc)
shclust <- hclust(as.dist(1 - ConsensusMatrix(stab)), method = stab$methods$linkage)
shclust <- as.hclust(rev(as.dendrogram(shclust)))

{
  pdf("Working_figures/Consensus_clustering_dendrogram.pdf",
    width = 15, height = 6
  )
  par(mar = c(1, 5, 1, 1))
  plot(as.dendrogram(shclust),
    xaxt = "n", leaflab = "none", axes = FALSE,
    xlab = "", sub = "", main = "", ylab = ""
  )
  axis(side = 2, at = seq(0, 1, by = 0.2), labels = NA)
  par(xpd = TRUE)
  text(-length(shclust$order) * 0.05,
    seq(0, 1, by = 0.2),
    seq(0, 1, by = 0.2),
    srt = 90
  )
  text(-length(shclust$order) * 0.08,
    mean(range(shclust$height)),
    "1 - co-membership proportion",
    srt = 90, cex = 1.5
  )
  rect.hclust(shclust, k = hat_N, border = "darkred")
  dev.off()
}

# Comparison
A <- matrix(0, nrow = length(myhclust$order), ncol = length(myhclust$order))
rownames(A) <- myhclust$order
colnames(A) <- shclust$order
for (i in 1:nrow(A)) {
  hm <- cutree(myhclust, k = hat_N)
  hm <- names(hm[which(hm == hm[paste0("obs", i)])])

  shm <- cutree(shclust, k = hat_N)
  shm <- names(shm[which(shm == shm[paste0("obs", i)])])
  overlap <- intersect(hm, shm)
  if (length(overlap) / length(union(hm, shm)) > 0.5) {
    A[as.character(i), as.character(i)] <- 1
  } else {
    A[as.character(i), as.character(i)] <- 2
  }
}
rownames(A) <- paste0("h", rownames(A))
colnames(A) <- paste0("sh", colnames(A))
Abig <- Square(A)

mycolours <- lighten(c("skyblue", "red", "forestgreen"), amount = 0.5)
myedgecolours <- darken(c("lightgrey", "red"), amount = 0.1)
g <- Graph(Abig, weighted = TRUE)
V(g)$label <- gsub("h", "", gsub("sh", "", V(g)$label))
V(g)$color <- mycolours[simul$theta[as.numeric(V(g)$label)]]
V(g)$size <- 6
V(g)$frame.color <- V(g)$color
V(g)$label.color <- darken(V(g)$color, amount = 0.5)
E(g)$color <- myedgecolours[E(g)$width]
mylayout <- matrix(c(
  rep(0, length(myhclust$order)), rep(1, length(myhclust$order)),
  1:length(myhclust$order), rev(1:length(myhclust$order))
), ncol = 2)

myasp <- 0.3
{
  pdf("Working_figures/Bipartite_network_comparison.pdf",
    width = myasp * 11, height = 10
  )
  par(mar = rep(0, 4))
  plot(g, layout = mylayout, asp = 1 / myasp)
  dev.off()
}

ClusteringPerformance(theta = cutree(myhclust, k = 3), theta_star = simul)
ClusteringPerformance(theta = cutree(shclust, k = hat_N), theta_star = simul)
