rm(list = ls())

library(sharp)

# Loading the results
stab_lc <- readRDS(paste0("Results/Inflammatory_proteins/Multivariate_lung_cancer.rds"))
expl_lc <- readRDS(paste0("Results/Inflammatory_proteins/Incremental_performances_lung_cancer.rds"))
stab_adeno <- readRDS(paste0("Results/Inflammatory_proteins/Multivariate_adenocarcinoma.rds"))
expl_adeno <- readRDS(paste0("Results/Inflammatory_proteins/Incremental_performances_adenocarcinoma.rds"))
stab_sclc <- readRDS(paste0("Results/Inflammatory_proteins/Multivariate_small-cell_carcinoma.rds"))
expl_sclc <- readRDS(paste0("Results/Inflammatory_proteins/Incremental_performances_small-cell_carcinoma.rds"))
smoking <- readRDS("Results/Inflammatory_proteins/Univariate_packyears.rds")

for (subtype in c("lc", "adeno", "sclc")) {
  tmpstab <- eval(parse(text = paste0("stab_", subtype)))
  tmpexpl <- eval(parse(text = paste0("expl_", subtype)))
  ids_selected <- tmpexpl$names[2:(sum(SelectedVariables(tmpstab)) + 1)]
  tmpsign <- apply(tmpexpl$Beta[[sum(SelectedVariables(tmpstab)) + 2]], 2,
    FUN = function(x) {
      prop.table(table(factor(sign(x), levels = c("-1", "1"))))
    }
  )[2, -c(1:3), drop = FALSE]
  tmpsign <- tmpsign[1, ids_selected]
  selprop <- SelectionProportions(tmpstab)[ids_selected]
  tmptable <- cbind(tmpsign, selprop)
  assign(paste0("mytable_", subtype), tmptable)
}

mytable <- cbind(rownames(mytable_lc),
  mytable_lc,
  data.frame(mytable_adeno)[rownames(mytable_lc), ],
  data.frame(mytable_sclc)[rownames(mytable_lc), ],
  beta = formatC(smoking[rownames(mytable_lc), "coef"], format = "f", digits = 2),
  pval = formatC(smoking[rownames(mytable_lc), "pval"], format = "e", digits = 2)
)
mytable <- as.matrix(mytable)
mytable[which(is.na(mytable))] <- ""
write.table(mytable, paste0("Tables/Inflammatory_proteins/Stably_selected.txt"),
  row.names = FALSE, col.names = TRUE, quote = FALSE, eol = "££\n", sep = "&"
)
