library(qtl)
library(MASS)
library(MESS)
library(ggplot2)
library(knitr)
source("haleyknott.R") # Linear mixed model with Haley-Knott regression 
source("utils.R")

ensure_directory("figures")
ensure_directory("logs")
log <- make_logger("logs/qtl_analysis_notes.md")

# load data 
WB02 <- read.cross(file="derived_data/Rqtl_CC27xC3H_BC.csv", na.strings=c("-", "NA", 'na'), 
                   genotypes=c("AA", "AB"), alleles=c("A", "B"), crosstype=c("bc"))
WB02 <- jittermap(WB02)
log("Cross object summary:")
#log(toString(summary(WB02)))

numF2s <- nrow(WB02$pheno)

# Plot genetic map
png("figures/genetic_map.png")
plotMap(WB02)
dev.off()

# Plot histograms of temperature data
png("figures/temp_histograms.png")
par(mfrow=c(3,2))
for (i in 6:10){
  plotPheno(WB02, pheno.col=i, xlab=expression(paste("Temperature (", degree, "C)")))
}
dev.off()

# Plot temperature data as time-series 
png("figures/temp_time-series.png")
xdata <- c(0,15,30,45,60) # Time post-challenge
ydata <- cbind(WB02$pheno[,'Temp_0min'], WB02$pheno[,'Temp_15min'], 
               WB02$pheno[,'Temp_30min'], WB02$pheno[,'Temp_45min'],
               WB02$pheno[,'Temp_60min'])
colors <- sample(rainbow(numF2s))
plot(xdata, ydata[1,], type = "l", col = colors[1], lty = 1, xlab = "Time post-challenge", 
     ylab="Temperature", main="Hypothermia in CC27xC3H backcross", ylim = c(24,40))
for (i in 2:numF2s){ # Number of mice 
  lines(xdata, ydata[i,], col=colors[i])
}
dev.off()

# Calculate averages
T00.avg = mean(ydata[,1])
T15.avg = mean(ydata[,2])
T30.avg = mean(ydata[,3])
T45.avg = mean(ydata[,4])
T60.avg = mean(ydata[,5])
avgs = c(T00.avg, T15.avg, T30.avg, T45.avg, T60.avg)

# Calculate summary statistic - area under the curve (technically area above the curve) 
Temp_aac <- rep(0,nrow(ydata))
for (i in 1:nrow(ydata)){
  yline <- ydata[i,]
  auc <- auc(x=xdata, y=yline)
  aac <- avgs[1]*60 - auc
  Temp_aac[i] <- aac 
}
WB02$pheno[,'Temp_aac'] <- Temp_aac

png("figures/temp_aac_hist.png")
plotPheno(WB02, pheno.col = 'Temp_aac')
dev.off()

# Data transformation 
WB02$pheno[,'Temp_aac'] <- trint(WB02$pheno[,'Temp_aac'])

png("figures/temp_aac_transformed_hist.png")
plotPheno(WB02, pheno.col = "Temp_aac")
dev.off()

# QTL analysis 
model <- hk(WB02, pheno.col = "Temp_aac")

png("figures/temp_aac_genome_scan.png")
plot(model, alternate.chrid = TRUE, ylab = "LOD", main = "Temperature - area above the curve")
dev.off()





