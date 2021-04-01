library(ggfortify)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}

organisms <- c('sce','kma','yli')
for (i in 1:length(organisms)){
  organism <- organisms[i]
  if (i==2){conditions <- c('Std','HiT','LpH')}
  else{conditions <- c('Std','HiT','LpH','Osm')}
#Load data
dataTypes <- c('fluxDist','absUsage')
for (j in 1:length(dataTypes)){
  dataType <- dataTypes[j]
  fileName <- paste('../../results/',organism,'_',dataType,'.txt',sep='')
  df  <- read.csv(file = fileName, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  #create numerical matrix of fluxes
  fluxDist <- as.data.frame(t(df[,5:ncol(df)]))
  if (dataType=='fluxDist'){colnames(fluxDist) <- df$rxns}
  if (dataType=='absUsage'){colnames(fluxDist) <- df$enzymes}
  #fluxDist[abs(fluxDist)==1E-15]<- 0
  fluxDist$info <- rownames(fluxDist)
  fluxDist <- separate(data = fluxDist, col = info, into = c("condition", "model"), sep = "_")
  fluxDist$condition <- factor(fluxDist$condition,levels=c('Std','HiT','LpH','Osm'))
  #fluxDist <- fluxDist[fluxDist$model!='GEM',]
  PCAdata  <- prcomp(fluxDist[,(1:(ncol(fluxDist)-2))], center = TRUE,scale = FALSE,retx=TRUE)
  p        <- autoplot(PCAdata,data = fluxDist,shape= 'model',colour = 'condition',size = 6,frame = TRUE, frame.type = 'norm')
  p        <- p + theme_bw(base_size = 24)
  p        <- p + scale_color_manual(values=c(rgb(0.45,0.45,0.45),rgb(0.7,0,0.3),rgb(0.8,0.5,0),rgb(0.1,0,0.9)))
  #p        <- p + stat_ellipse(aes(group=model),level=0.99)
  plotName <- paste('../../results/Figure_S2/',organism,'_',dataType,'_PCA.png',sep='') 
  png(plotName,width = 650, height = 600)
  plot(p)
  dev.off()
}
}
