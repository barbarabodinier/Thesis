rm(list=ls())
setwd("~/Dropbox/PhD/Thesis/Version2/chapter1/")

library(data.table)
library(rworldmap)

# Reading latest global age-standardised rates per cancer type
# Downloaded from https://gco.iarc.fr/today/online-analysis-multi-bars?v=2020&mode=cancer&mode_population=countries&population=900&populations=900&key=asr&sex=0&cancer=39&type=0&statistic=5&prevalence=0&population_group=0&ages_group%5B%5D=0&ages_group%5B%5D=17&nb_items=55&group_cancer=0&include_nmsc=1&include_nmsc_other=0&type_multiple=%257B%2522inc%2522%253Atrue%252C%2522mort%2522%253Atrue%252C%2522prev%2522%253Afalse%257D&orientation=horizontal&type_sort=0&type_nb_items=%257B%2522top%2522%253Atrue%252C%2522bottom%2522%253Afalse%257D
cancer_types=read.table("Data/asr_cancer_types_worldwide_globocan_2020.csv", sep=",")
cancer_types=cancer_types[,-ncol(cancer_types)]
colnames(cancer_types)=c("Incidence", "Mortality")

# Identifying sex-stratified rates
female_specific=c(1,4,9,11,30,35)
male_specific=c(2,24,31)
rownames(cancer_types)[female_specific]=paste0(rownames(cancer_types)[female_specific], "*")
rownames(cancer_types)[male_specific]=paste0(rownames(cancer_types)[male_specific], "**")

# Parameters
mylwd=7
mycolours=c("navy", "red")

plotname="Working_figure/Incidence_mortality_by_cancer_type_worldwide.pdf"
pdf(plotname, width = 12, height=7)
par(mar=c(12.5, 5, 1, 1))
plot(seq(1, nrow(cancer_types))+0.1,
     cancer_types[,1], 
     xlim=c(1.5, nrow(cancer_types)+1),
     ylim=c(0, 50),
     type="h", lend=1, lwd=mylwd,
     bty="n",
     xlab="", xaxt="n", 
     ylab="Age-standardised rate per 100,000",
     cex.lab=1.25,
     col=mycolours[1])
points(seq(1, nrow(cancer_types))+0.4,
       cancer_types[,2], 
       type="h", lend=1, lwd=mylwd,
       col=mycolours[2])
axis(side=1, at=seq(1, nrow(cancer_types))+0.25, 
     labels = rownames(cancer_types), las=2)
dev.off()
system(paste("pdfcrop --margin 10",plotname,plotname))

plotname="Working_figure/Incidence_mortality_by_cancer_type_worldwide_legend.pdf"
pdf(plotname, width = 4, height=4)
par(mar=rep(5, 4))
plot.new()
legend("top", pch=15, col=mycolours, legend=c("Incidence", "Mortality"), bty="n")
dev.off()
system(paste("pdfcrop --margin 10",plotname,plotname))

# Reading country-specific age-standardised mortality of lung cancer
# Downloaded from https://gco.iarc.fr/today/online-analysis-map?v=2020&mode=population&mode_population=continents&population=900&populations=900&key=asr&sex=0&cancer=15&type=1&statistic=5&prevalence=0&population_group=0&ages_group%5B%5D=0&ages_group%5B%5D=17&nb_items=10&group_cancer=1&include_nmsc=0&include_nmsc_other=0&projection=natural-earth&color_palette=default&map_scale=quantile&map_nb_colors=5&continent=5&show_ranking=0&rotate=%255B10%252C0%255D
spatial_lung=read.table("Data/mortality_lung_cancer_europe_globocan_2020.csv", sep=",")
spatial_lung=cbind(rownames(spatial_lung), spatial_lung)
colnames(spatial_lung)=c(colnames(spatial_lung)[-1], "NA")
spatial_lung=spatial_lung[,-ncol(spatial_lung)]

pdf("Working_figure/Map_europe_lung_cancer_mortality.pdf",
    width = 10, height=7)
par(mar=rep(5,4))
spdf=joinCountryData2Map(dF=spatial_lung, 
                         joinCode = "ISO3", 
                         nameJoinColumn="ISO.code", 
                         mapResolution = "low")
mapCountryData(spdf, 
               mapTitle = "",
               lwd=1,
               catMethod = "pretty",
               numCats = 100,
               addLegend = TRUE,
               colourPalette = colorRampPalette(c("white", "darkred"))(100),
               borderCol = "white",
               nameColumnToPlot = "Value",
               missingCountryCol = "grey",
               xlim = c(15, 20),
               ylim = c(37, 70))
dev.off()

# Reading age-standardised rates of lung cancer over time in Norway and Italy
# Downloaded from https://gco.iarc.fr/overtime/en/dataviz/trends?populations=57800_38000&sexes=1_2&types=1&multiple_cancers=0&mode=population&multiple_populations=1&smoothing=0&scale=linear&min_zero=1
temporal_lung=fread("Data/mortality_trends_over_time_lung_cancer_globocan_2020.csv", data.table = FALSE)
rownames(temporal_lung)=temporal_lung[,1]
temporal_lung=temporal_lung[,-1]

mypch=c(17, 19, 17, 19)
mycolours=c(rep("forestgreen", 2), rep("chocolate", 2))

pdf("Working_figure/Time_trends_lung_cancer_mortality.pdf",
    width = 10, height=7)
par(mar=rep(5,4))
plot(NA, 
     xlim=range(as.numeric(temporal_lung[1,]), na.rm = TRUE), 
     ylim=c(0,60), las=1,
     xlab="Year", cex.lab=1.25,
     ylab="Age-standardised mortality rate per 100,000",
     panel.first=abline(v=seq(1950, 2020,by=10), lty=3, col="grey"))
for (k in 2:5){
  points(as.numeric(temporal_lung[1,]), 
         as.numeric(temporal_lung[k,]), 
         pch=mypch[k-1], 
         cex=0.5,
         col=mycolours[k-1])
  lines(as.numeric(temporal_lung[1,]), 
        as.numeric(temporal_lung[k,]),
        col=mycolours[k-1])
}
dev.off()

