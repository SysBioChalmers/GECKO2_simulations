library(dplyr)
library(fmsb)

if (exists("RStudio.Version")){
  current <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(current)
} else {
  setwd(getSrcDirectory()[1])
} 
#Load data
fileName <- '../../results/Ecoli_Ptot_comparison.txt'
expData  <- read.csv(file = fileName, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

#plot
maxLim  <- maxLim
minLim  <- 0
#generateSpiderPlots(exchData,feature,'x',0,1)
plotName    <- '../../results/Figure_2/Ecoli_Ptot_predictions.png'

#generateSpiderPlots <- function(df,feature,exclude,minLim,maxLim){
maxLim <- 0.8
minLim <- 0.4
plotDf <- (t(expData[,2:ncol(expData)]))
#rownames(plotDf) <- c('')
plotDf <- rbind(rep(maxLim,ncol(plotDf)),rep(minLim,ncol(plotDf)),plotDf)
plotDf <- as.data.frame(plotDf)
colnames(plotDf) <- expData$cSource
#==================
# Plot 2: Same plot with custom features
colors_border = c(rgb(0.8,0.6,0,0.8), rgb(0.1,0,0.8,0.8), rgb(0.4,0.4,0.40,0.8))
colors_in     = c(rgb(0.8,0.6,0,0.2), rgb(0.1,0,0.8,0.2), rgb(0.4,0.4,0.40,0.2))
png(plotName,width = 650, height = 500)
radarchart( plotDf , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(minLim,maxLim*100,(maxLim-minLim)/4), cglwd=1.5,
            #custom labels
            vlcex=1.5, calcex = 1.5)
#plot(p)
legend(x=1.05, y=1.20, legend = colnames(expData)[2:ncol(expData)], bty = "n", pch=20 , col=colors_in , text.col = "black", cex=1.5, pt.cex=3)
dev.off()

#}




