library(reshape)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
#Get repo directory
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
organisms <- c('sce','kma','yli')
for (organism in organisms){
  print(organism)
  conditions <- c('Std','HiT','LpH','Osm')
  if (all(organism=='yli')){conditions <- c('Std','HiT','LpH')}
  models     <- c('GEM','ecM','ecP')
  repoPath <- getwd()
  repoPath <- gsub('/code/analyse_protModels','',repoPath)
  setwd(repoPath)
  fileName <- paste(repoPath,'/results/',organism,'_fluxDist.txt',sep='')
  dataset  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
 #Discard those fluxes with a zero value for all conditions
  dataset  <- dataset[rowSums(dataset[,5:ncol(dataset)])>0,]
    
  data_GEM <- dataset[,c(1,grep('GEM',colnames(dataset)))]
  data_GEM <- melt(data_GEM, id='rxns')
  data_GEM <- separate(data_GEM,'variable',c('condition','model'),sep='_')
  
  data_ecM  <- dataset[,c(1,grep('ecM',colnames(dataset)))]
  data_ecM  <- melt(data_ecM, id='rxns')
  data_ecM  <- separate(data_ecM,'variable',c('condition','model'),sep='_')
  
  data_ecP  <- dataset[,c(1,grep('ecP',colnames(dataset)))]
  data_ecP  <- melt(data_ecP, id='rxns')
  data_ecP  <- separate(data_ecP,'variable',c('condition','model'),sep='_')
  
  df1 <- data.frame(data_GEM$rxns,data_GEM$condition,as.numeric(data_GEM$value)+1E-7)
  df1 <- rbind(df1,df1)
  
  x1 <- as.data.frame(data_ecM$value+1E-7)
  x1 <- cbind(x1,rep('ecM',nrow(x1)))
  colnames(x1) <-c('ecModels','model')
  
  x2 <- as.data.frame(data_ecP$value+1E-7)
  x2 <- cbind(x2,rep('ecP',nrow(x2)))
  colnames(x2) <-c('ecModels','model')
  
  df2 <- rbind(x1,x2)
  df <- c()
  df  <- cbind(df1,df2)
  df  <- cbind(df1,df2)
  colnames(df) <- c('rxns','condition','GEM','ecModels','models')
  #df <- df[rowSums((df[,3:4]))>1.9E-7,]
  df$FC <- log10((df$ecModels)/(df$GEM+1E-7))
  #discard those fluxes that change direction
  df <- df[!is.na(df$FC),]
  df$Fchange <- rep('good',nrow(df))
  df$Fchange[abs(df$FC)>log10(2)] <- 'bad'
  df$Fchange[abs(df$FC)>1] <- 'worse'
  df$Fchange <- (df$Fchange)
  df$condition <- factor(df$condition,levels=conditions)
  #create table with results summary
  conditions <- factor(unique(df$condition),levels = c('Std','HiT','LpH','Osm'))
  models     <- unique(df$models)
  predLvls   <- c('good','bad','worse')
  newTable <- data.frame()
  for (cond in conditions){
    total <- c()
    for (mType in models){
      total  <- intersect(grep(cond,df$condition),grep(mType,df$models))
      #get fluxes that are on/off for the ecModels
      subDF <- df[total,]
      flux_ON  <- which(subDF$GEM==1E-7 & subDF$ecModels>1E-7 & subDF$models == mType)
      flux_OFF <- which(subDF$GEM>1E-7 & subDF$ecModels==1E-7 & subDF$models == mType)
      newRow <- data.frame(cond,mType,length(total),stringsAsFactors = TRUE)
      for (predLevel in predLvls){
        positions <- intersect(total,grep(predLevel,df$Fchange))
        fraction  <- length(positions)/length(total)
        newRow    <- cbind(newRow,fraction)
      }
      fractionON  <- length(flux_ON)/length(total)
      fractionOFF <- length(flux_OFF)/length(total)
      newRow      <- cbind(newRow,fractionON,fractionOFF)
      newTable    <- rbind(newTable,newRow)
    }
  }
  colnames(newTable) <- c('condition','model','total',predLvls[1],predLvls[2],predLvls[3],'fluxes_ON','fluxes_OFF')
  
  fileName <- paste(repoPath,'/results/',organism,'_fluxComp_summary.txt',sep='')
  write.table(newTable, file = fileName, row.names = FALSE,quote= FALSE,sep = '\t')
  colourCount = length(conditions)
  getPalette = colorRampPalette(brewer.pal(5, "Set1"))
  
  p <- ggplot(df, aes(x=GEM, y=ecModels, shape=models, color=condition)) +
    geom_point(size=4) +
    scale_x_log10(limits=c(1E-7,1E2)) + scale_y_log10(limits = c(1E-7, 1E2)) +
    theme_bw(base_size = 32) + scale_color_manual(values=c(rgb(0.45,0.45,0.45),rgb(0.7,0,0.3),rgb(0.8,0.6,0),rgb(0.1,0,0.9))) +
    geom_line(data=df,aes(x=GEM,y=GEM*2),linetype="dashed",color = "dark grey", size=0.5) +
    geom_line(data=df,aes(x=GEM,y=GEM*0.5),linetype="dashed",color = "dark grey", size=0.5)+
    geom_line(data=df,aes(x=GEM,y=GEM*10),linetype="dashed",color = "black", size=0.5) +
    geom_line(data=df,aes(x=GEM,y=GEM*0.1),linetype="dashed",color = "black", size=0.5) +
    xlab('Predicted fluxes (GEM) [mmol/gDW h]') + ylab('Predicted fluxes (ecModels) [mmol/gDw h]') +
    labs(shape="ecModels", colour="Conditions")
  plotName <- paste(repoPath,'/results/Figure_S2/',organism,'_pairwiseFluxComp.png',sep='')
  print(plotName)
  png(plotName,width = 850, height = 750)
  plot(p)
  dev.off()
}