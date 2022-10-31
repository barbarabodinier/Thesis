rm(list = ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter4/")

library(fake)
library(sharp)
library(igraph)

source("Scripts/additional_functions_specific_to_comparisons_stability.R")

# Simulation parameters
simul_study_id <- 1
topology <- "random"
PFER_thr <- 20

# Saving figure
metric_list <- c("precision", "recall", "F1_score", "time.user.self")
mytable=NULL
for (simul_id in 1:3) {
  performances <- readRDS(paste0("Results/Stability_selection/2-simulations/1-graphical_model/Sensitivity_K_", simul_study_id, "_", topology, "/Performances_", simul_id, "_merged.rds"))
  performances <- performances[, , 1:min(dim(performances)[3], 1000)]

  tmptable <- matrix(paste0(
    formatC(apply(performances[, metric_list, ], c(1, 2), median), format = "f", digits = 3),
    " [",
    formatC(apply(performances[, metric_list, ], c(1, 2), IQR), format = "f", digits = 3),
    "]"
  ),
  ncol = length(metric_list)
  )
  tmptable=cbind(rep("", nrow(tmptable)),
                 c(10, 20, 50, 100, 500, 1000, 2000, 5000), 
                 tmptable)
  colnames(tmptable)=c("Dimensionality","K", metric_list)
  
  mytable=rbind(mytable, tmptable)
}
write.table(mytable, paste0("Tables/Table_precision_recall_K_", simul_study_id, "_", topology, ".txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, eol = "££\n", sep = "&"
)
