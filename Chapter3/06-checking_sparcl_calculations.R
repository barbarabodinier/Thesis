rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter3/")

library(sharp)
library(sparcl)

# Generate 2-class data
set.seed(1)
x <- matrix(rnorm(100*50),ncol=50)
y <- c(rep(1,50),rep(2,50))
x[y==1,1:25] <- x[y==1,1:25]+2
# Do tuning parameter selection for sparse hierarchical clustering
perm.out <- HierarchicalSparseCluster.permute(x, wbounds=c(1.5,2:6),
                                              nperms=5)
print(perm.out)
plot(perm.out)
# Perform sparse hierarchical clustering
sparsehc <- HierarchicalSparseCluster(dists=perm.out$dists,
                                      wbound=perm.out$bestw, method="complete")
as.matrix(dist(x))[1:5,1:5]
sparsehc$u[1:5,1:5]
sparsehc$ws

mydists=list()
mydist=matrix(0, nrow=nrow(x), ncol=nrow(x))
for (k in 1:ncol(x)){
  mydists=c(mydists, list(as.matrix(dist(x[,k,drop=FALSE]))))
  tmpdist=as.matrix(dist(x[,k,drop=FALSE]))^2
  mydist=mydist+sparsehc$ws[k]*tmpdist
}
mydist[1:5,1:5]

plot(sparsehc$u[upper.tri(sparsehc$u)], mydist[upper.tri(mydist)])


