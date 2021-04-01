library(dplyr)
library(fmsb)

if (exists("RStudio.Version")){
  current <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(current)
} else {
  setwd(getSrcDirectory()[1])
}

df         <- data.frame(stringsAsFactors = FALSE)
organisms  <- c('sce','kma','yli')
orgKeys    <- c('ecYeastGEM','eciSM966','eciYali')

for (index in 1:length(organisms)){
  organism <- organisms[index]
  orgKey   <- orgKeys[index]
  conditions <- c('Std','HiT','LpH','Osm')
  maxLim <- 1
  if (all(organism=='yli')){
    conditions <- c('Std','HiT','LpH')
    maxLim <- 0.3
  }
  #Load data
  fileName <- paste('../',organism,'_scripts/fermentationData.txt',sep='')
  expData  <- read.csv(file = fileName, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  dataPath <- paste('../../ecModels/',orgKey,'_prot/',sep='')
  #open exchange fluxes data for GEM, ecModel and ecModel_prot
  fileName    <- paste('../../results/',organism,'_exchFluxes_GEM.txt',sep='')
  fluxes_gem  <- read.csv(file = fileName, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  fileName    <- paste('../../results/',organism,'_exchFluxes_ecM.txt',sep='')
  fluxes_ecM  <- read.csv(file = fileName, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  fileName    <- paste('../../results/',organism,'_exchFluxes_pro.txt',sep='')
  fluxes_pro  <- read.csv(file = fileName, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  error_df    <- c()#data.frame(stringsAsFactors = FALSE)
  for (j in 1:length(conditions)){
    cond <- conditions[j]
    error <- 0
    fluxesExp <- c(expData$Drate[j],expData$oxygen[j],expData$CO2[j])
    fluxes    <- c(fluxes_gem$Drate[j],fluxes_gem$oxygen[j],fluxes_gem$CO2[j])
    error_gem <- mean(abs(abs(fluxesExp)-abs(fluxes))/abs(fluxesExp))
    fluxes    <- c(fluxes_ecM$Drate[j],fluxes_ecM$oxygen[j],fluxes_ecM$CO2[j])
    error_ecM <- mean(abs(abs(fluxesExp)-abs(fluxes))/abs(fluxesExp))
    fluxes    <- c(fluxes_pro$Drate[j],fluxes_pro$oxygen[j],fluxes_pro$CO2[j])
    error_pro <- mean(abs(abs(fluxesExp)-abs(fluxes))/abs(fluxesExp))
    errors    <- c(error_pro,error_gem,error_ecM)
    error_df  <- cbind(error_df,errors)
  }
  rownames(error_df) <- c('ecP','GEM','ecM')  
  colnames(error_df) <- conditions
  error_df <- as.data.frame(error_df,stringsAsFactors = FALSE)
  #plot
  maxLim  <- maxLim
  minLim  <- 0
  #generateSpiderPlots(exchData,feature,'x',0,1)
  plotName    <- paste('../../results/Figure_4/exchFluxes_',organism,'.png',sep='')
  
  #generateSpiderPlots <- function(df,feature,exclude,minLim,maxLim){
  plotDf <- (error_df)
  plotDf <- rbind(rep(maxLim,ncol(plotDf)),rep(minLim,ncol(plotDf)),plotDf)
  #==================
  # Plot 2: Same plot with custom features
  colors_border = c(rgb(0.8,0.6,0,0.8), rgb(0.1,0,0.8,0.8), rgb(0.4,0.4,0.40,0.8))
  colors_in     = c(rgb(0.8,0.6,0,0.2), rgb(0.1,0,0.8,0.2), rgb(0.4,0.4,0.40,0.2))
  png(plotName,width = 550, height = 500)
  radarchart( plotDf , axistype=1 , 
              #custom polygon
              pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(minLim,maxLim*100,(maxLim-minLim)/4), cglwd=1.5,
              #custom labels
              vlcex=2, calcex = 1.5)
  #plot(p)
  legend(x=0.85, y=0.85, legend = rownames(error_df), bty = "n", pch=20 , col=colors_in , text.col = "black", cex=1.5, pt.cex=3)
  dev.off()
}
  #}
  
  
