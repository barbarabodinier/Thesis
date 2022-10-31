rm(list = ls())

library(fake)

# Simulation of data with clusters
set.seed(0)
n <- c(5, 15, 10)
simul <- SimulateClustering(
  n = n,
  pk = 20,
  nu_xc = 1,
  ev_xc = 0.3
)
x <- simul$data

# Hierarchical clustering
mydist <- as.matrix(dist(x))
myhclust <- hclust(dist(x), method = "complete")
hclustorder <- myhclust$order

pdf(paste0("Working_figures/Hierarchical_clustering_example.pdf"),
  width = 10, height = 7
)
mycolours <- c("royalblue", "tomato", "tan")
par(mar = c(5, 5, 1, 1))
plot(as.dendrogram(myhclust),
  xaxt = "n", leaflab = "none", axes = FALSE,
  xlab = "", sub = "", main = "",
  ylab = "Euclidian distance", cex.lab = 1.5
)
axis(side = 2, at = seq(0, max(myhclust$height), by = 2), las = 1)
rect.hclust(myhclust, k = 3, border = "navy")
for (k in 1:length(myhclust$order)) {
  axis(
    side = 1, at = k, paste0("X", myhclust$order[k]),
    las = 2,
    col.axis = mycolours[cutree(myhclust, k = 3)[myhclust$order[k]]]
  )
}
dev.off()
