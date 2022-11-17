rm(list = ls())

library(openxlsx)
library(sharp)
library(summarise)

# Loading the results
stab <- readRDS("Results/2-clustering/consensus_clustering_very_fine_both.rds")

# Extracting the clusters
myclusters <- Clusters(stab)

# Extracting the multivariate clusters (>1 members)
multiv_clusters <- sort(unique(myclusters[duplicated(myclusters)]))
mymultiv_clusters <- myclusters[myclusters %in% multiv_clusters]

# Extracting mass and retention time
mass <- as.numeric(gsub("@.*", "", names(mymultiv_clusters)))
time <- as.numeric(gsub(".*@", "", names(mymultiv_clusters)))

# Calculating the range by multivariate cluster
for (characteristic in c("time", "mass")) {
  print(characteristic)
  for (method in c("min_max", "inter_quartile")) {

    # Calculating the difference between max and min value per cluster
    mylist <- split(get(characteristic), f = mymultiv_clusters)
    mylist <- mylist[sort.list(sapply(mylist, length))]
    if (method == "min_max") {
      myranges <- sapply(mylist, FUN = function(x) {
        diff(range(x))
      })
    } else {
      myranges <- sapply(mylist, FUN = function(x) {
        IQR(x)
      })
    }

    # Preparing figure
    pdf(paste0("Figures/2-clustering/Range_", characteristic, "_", method, ".pdf"))
    par(mar = c(5, 5, 1, 1))
    tmp <- hist(myranges,
      breaks = seq(
        from = 0,
        to = ifelse(characteristic == "time", yes = 7, no = 1000),
        by = ifelse(characteristic == "time", yes = 0.1, no = 10)
      ),
      xlab = paste0(
        ifelse(method == "min_max",
          yes = "Range in ",
          no = "Inter-quartile range in "
        ),
        ifelse(characteristic == "time",
          yes = "retention time (s)",
          no = "mass"
        )
      ),
      ylab = "Counts", cex.lab = 1.5, main = "",
      col = "navy", border = "white"
    )
    dev.off()

    # Preparing table
    mybreaks <- paste0("]", tmp$breaks[-length(tmp$breaks)], ", ", tmp$breaks[-1], "]")
    mybreaks[1] <- gsub("]", "[", mybreaks[1])
    mytable <- cbind(Range = mybreaks, Count = tmp$counts)
    write.xlsx(as.data.frame(mytable),
      paste0("Tables/2-clustering/Range_", characteristic, "_", method, ".xlsx"),
      rowNames = FALSE, colNames = TRUE
    )
  }
}
