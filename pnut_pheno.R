### Phenotype processing functions for peanut allergy data 
### phenotypes - temperature, symptom score 

#----------------------------Phenotype plotting--------------------------------#

# Temperature time-series 
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

# Averages 
T00.avg = mean(ydata[,1])
T15.avg = mean(ydata[,2])
T30.avg = mean(ydata[,3])
T45.avg = mean(ydata[,4])
T60.avg = mean(ydata[,5])
avgs = c(T00.avg, T15.avg, T30.avg, T45.avg, T60.avg)
points(xdata, avgs, type="l", col="black")

#-----------------------------Derived measures---------------------------------#

Temp_aac <- rep(0,nrow(ydata))

for (i in 1:nrow(ydata)){
  yline <- ydata[i,]
  auc <- auc(x=xdata, y=yline)
  aac <- avgs[1]*60 - auc
  Temp_aac[i] <- aac 
}

WB02$pheno[,'Temp_aac'] <- Temp_aac