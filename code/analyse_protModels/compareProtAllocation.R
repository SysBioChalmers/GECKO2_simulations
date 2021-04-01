library(reshape)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(matrixStats)
library(pheatmap)
library(viridis)
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
  fileName <- paste(repoPath,'/results/',organism,'_absUsage.txt',sep='')
  dataset  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  #
  fileName <- paste(repoPath,'/results/',organism,'_relUsage.txt',sep='')
  dataset_rel  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  #Discard those fluxes with a zero value for all conditions
  nonZeros    <- rowSums(dataset[,5:ncol(dataset)])>0
  dataset     <- dataset[nonZeros,]
  dataset_rel <- dataset_rel[nonZeros,]
  #Reshape dataset 
  data_ecM <- dataset[,c(1,grep('ecM',colnames(dataset)))]
  data_ecM <- melt(data_ecM, id='enzymes')
  data_ecM <- separate(data_ecM,'variable',c('condition','model'),sep='_')
  data_ecP <- dataset[,c(1,grep('ecP',colnames(dataset)))]
  data_ecP <- melt(data_ecP, id='enzymes')
  data_ecP <- separate(data_ecP,'variable',c('condition','model'),sep='_')
  data_rel <- dataset_rel[,c(1,grep('ecP',colnames(dataset_rel)))]
  data_rel <- melt(data_rel, id='enzymes')
  data_rel <- separate(data_rel,'variable',c('condition','model'),sep='_')
  df <- data.frame(data_ecM$enzymes,data_ecM$condition,dataset$subSystems,as.numeric(data_ecM$value)+1E-12,as.numeric(data_ecP$value)+1E-12,as.numeric(data_rel$value))
  colnames(df) <- c('enzymes','condition','subSystems','ecM','ecP','rel')
  relDF <- df[df$rel!=Inf & !is.na(df$rel),]
  relDF$protAb <- relDF$ecP/relDF$rel
  relDF$protAb[relDF$protAb==Inf] <- NaN
  df$FC <- log10((df$ecP)/(df$ecM))
  #discard those fluxes that change direction
  df <- df[!is.na(df$FC),]
  df$precision <- rep('good',nrow(df))
  df$precision[abs(df$FC)>log10(2)] <- 'bad'
  df$precision[abs(df$FC)>1] <- 'worse'
  df$precision <- (df$precision)
  df$condition <- factor(df$condition,levels=conditions)
  relDF$condition <- factor(relDF$condition,levels=conditions)
  
  #create table with results summary
  conditions <- factor(unique(df$condition),levels = c('Std','HiT','LpH','Osm'))
  predLvls   <- c('good','bad','worse')
  newTable <- data.frame()
  #Get tables with clasification of predictions precision (good, bad, worse ecM vs ecP)
  for (cond in conditions){
    total <- c()
    total <- (grep(cond,df$condition))
    #get fluxes that are on/off for the ecModels
    subDF    <- df[total,]
    flux_ON  <- which(subDF$ecM==1E-12 & subDF$ecP>1E-12)
    flux_OFF <- which(subDF$ecM>1E-12 & subDF$ecP==1E-12)
    newRow <- data.frame(cond,length(total),stringsAsFactors = TRUE)
    for (predLevel in predLvls){
      positions <- intersect(total,grep(predLevel,df$precision))
      fraction  <- length(positions)/length(total)
      newRow    <- cbind(newRow,fraction)
    }
    fractionON  <- length(flux_ON)/length(total)
    fractionOFF <- length(flux_OFF)/length(total)
    newRow      <- cbind(newRow,fractionON,fractionOFF)
    newTable    <- rbind(newTable,newRow)
  }
  colnames(newTable) <- c('condition','total',predLvls[1],predLvls[2],predLvls[3],'fluxes_ON','fluxes_OFF')
  fileName <- paste(repoPath,'/results/',organism,'_enzUsageComp_summary.txt',sep='')
  write.table(newTable, file = fileName, row.names = FALSE,quote= FALSE,sep = '\t')
  colourCount = length(conditions)
  getPalette = colorRampPalette(brewer.pal(5, "Set1"))
  #Plot pairwise comparison of enzyme usage profiles across conditions
  p <- ggplot(df, aes(x=ecM, y=ecP,color=condition,shape=condition)) +
    geom_point(size=3) +
    scale_x_log10(limits=c(1E-12,1E-3)) + scale_y_log10(limits = c(1E-12, 1E-3)) +
    theme_bw(base_size = 26) + scale_color_manual(values=c(rgb(0.45,0.45,0.45),rgb(0.7,0,0.3),rgb(0.8,0.6,0),rgb(0.1,0,0.9))) +
    geom_line(data=df,aes(x=ecM,y=ecM*2),linetype="dashed",color = "dark grey", size=0.5) +
    geom_line(data=df,aes(x=ecM,y=ecM*0.5),linetype="dashed",color = "dark grey", size=0.5)+
    geom_line(data=df,aes(x=ecM,y=ecM*10),linetype="dashed",color = "black", size=0.5) +
    geom_line(data=df,aes(x=ecM,y=ecM*0.1),linetype="dashed",color = "black", size=0.5) +
    xlab('Enz usage (ecM) [mmol/gDW]') + ylab('Enz usage (ecP) [mmol/gDw]') +
    labs(colour="Conditions")
  plotName <- paste(repoPath,'/results/Figure_S2/',organism,'_pairwise_EnzUseComp.png',sep='')
  png(plotName,width = 650, height = 600)
  plot(p)
  dev.off()
  
  #Get tables with characterization of relative-absolute usage relations
  relDF$rel_int <- rep('medium',nrow(relDF))
  relDF$rel_int[abs(relDF$rel)<0.25] <- 'low'
  relDF$rel_int[abs(relDF$rel)>0.75] <- 'high'
  relDF$abs_int <- rep('medium',nrow(relDF))
  lowThrsld     <- 0.5*median(relDF$protAb,na.rm = TRUE)#-sd(relDF$protAb,na.rm = TRUE)
  highThrsld    <- 2*median(relDF$protAb,na.rm = TRUE)#+sd(relDF$protAb,na.rm = TRUE)
  relDF$abs_int[abs(relDF$protAb)<lowThrsld]  <- 'low'
  relDF$abs_int[abs(relDF$protAb)>highThrsld] <- 'high'
  relDF$condition <- factor(relDF$condition,levels=conditions)
  usgLvls         <- c('low','medium','high')
  relDF$abs_int   <- factor(relDF$abs_int,levels=usgLvls)
  relDF$rel_int   <- factor(relDF$rel_int,levels=usgLvls)
  newTable        <- data.frame(stringsAsFactors = TRUE)
  for (cond in conditions){
    #get fluxes that are on/off for the ecModels
    subDF <- relDF[grep(cond,relDF$condition),]
    total <- nrow(subDF)
    for (absLevel in usgLvls){
      newRow <- data.frame(cond,total,absLevel,stringsAsFactors = TRUE)
      for (relLevel in usgLvls){
        positions <- intersect(grep(absLevel,subDF$abs_int),grep(relLevel,subDF$rel_int))
        fraction  <- length(positions)/total
        newRow    <- cbind(newRow,fraction)
      }
      newTable <- rbind(newTable,newRow)
    }
  }
  colnames(newTable) <- c('condition','total','abs_usage','relUsg_low','relUsg_medium','relUsg_high')
  fileName <- paste(repoPath,'/results/',organism,'_usage_relVsAbs_summary.txt',sep='')
  write.table(newTable, file = fileName, row.names = FALSE,quote= FALSE,sep = '\t')
  
  # colourCount = length(conditions)
  # getPalette = colorRampPalette(brewer.pal(5, "Set1"))
  # #Plot pairwise comparison of absolute and relative protein usage for ecP
  # p <- ggplot(relDF, aes(x=protAb, y=rel,color=condition,shape=condition)) +
  #   geom_point(size=4) +
  #   scale_x_log10(limits=c(1E-7,1E-3)) + scale_y_continuous(limits = c(0, 1)) +
  #   theme_bw(base_size = 32) + scale_color_manual(values=c(rgb(0.45,0.45,0.45),rgb(0.7,0,0.3),rgb(0.8,0.6,0),rgb(0.1,0,0.9))) +
  #   geom_hline(data=relDF,aes(yintercept=0.25),linetype="dashed",color = "black", size=0.5) +
  #   geom_hline(data=relDF,aes(yintercept=0.75),linetype="dashed",color = "black", size=0.5)+
  #   geom_vline(data=relDF,aes(xintercept=highThrsld),linetype="dashed",color = "black", size=0.5) +
  #   geom_vline(data=relDF,aes(xintercept=lowThrsld),linetype="dashed",color = "black", size=0.5) +
  #   xlab('Protein abundances [mmol/gDW]') + ylab('Relative enzyme usages') +
  #   labs(colour="Conditions")
  # plotName <- paste(repoPath,'/results/Figure3_S2_S3/',organism,'_enzUsage_relVsAbs_Comp.png',sep='')
  # png(plotName,width = 850, height = 750)
  # plot(p)
  # dev.off()
  
  #Get saturated enzymes table
  fileName <- paste(repoPath,'/data/',organism,'_enzymeTable.txt',sep='')
  enzInfo  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  fileName <- paste(repoPath,'/data/model_enz_OGs.txt',sep='')
  modelOGs <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  fileName <- paste(repoPath,'/results/',organism,'_relUsage.txt',sep='')
  temp     <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  temp     <- temp[,c(1,2,3,4,grep('ecP',colnames(temp)))]
  satDF    <- temp
  #Extract orthologs
  satDF$OG <- enzInfo$OG[match(satDF$enzymes,enzInfo$enzymes)]
  satDF    <- satDF[!is.na(satDF$OG),]
  idxs     <- match(modelOGs$Ogs_models,satDF$OG)
  satDF    <- satDF[idxs[!is.na(idxs)],]
  fileName <- paste(repoPath,'/results/',organism,'_ALLOG_relUsage.txt',sep='')
  write.table(satDF, file = fileName, row.names = FALSE,quote= FALSE,sep = '\t')
  #Get stress related proteins (just keep the proteins that are used/measured in at least one condition)
  temp  <- do.call(data.frame,lapply(temp, function(x) replace(x, is.infinite(x),0)))
  idxs  <- which(rowSums(temp[,5:ncol(temp)])>0)
  temp  <- temp[idxs,]
  #get saturated proteins for at least one condition
  idxs  <- which(rowSums(temp[,5:ncol(temp)] >= 0.95) >=1)
  temp  <- temp[idxs,]
  #Get proteomics abundance
  fileName  <- paste(repoPath,'/results/ecModels_',organism,'_protConcs.txt',sep='')
  protConcs <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  idxs <- match(temp$enzymes,protConcs$enzymes)
  protConcs  <- protConcs[idxs,]
  protConcs  <- protConcs[,6:ncol(protConcs)]
  protConcs  <- protConcs/protConcs[,1]
  matriz <- temp[,5:ncol(temp)]*0
  for (i in 1:nrow(matriz)){
    pos1 <-  which(protConcs[i,]>=0.95)
    pos2 <-  which(temp[i,5:ncol(temp)]>=0.95)
    matriz[i,intersect(pos1,pos2)] <- 1
  }
  toKeep <- rowSums(matriz)>0
  temp[,5:ncol(temp)] <- matriz
  temp <- temp[toKeep,]
  #Write file and get heatmap
  idxs <- match(temp$enzymes,enzInfo$enzymes)
  temp$OG <- enzInfo$OG[idxs]
  fileName <- paste(repoPath,'/results/',organism,'_sat_Enz.txt',sep='')
  write.table(temp, file = fileName, row.names = FALSE,quote= FALSE,sep = '\t')
  matriz   <- temp
  rownames(matriz) <- matriz$enzNames
  idxs <- which(matriz$Std_ecP==1 & rowSums(matriz[,5:(ncol(matriz)-1)])<length(conditions))
  matriz <- matriz[-idxs,]
  matriz <- matriz[,5:(ncol(matriz)-1)]
  colnames(matriz) <- conditions
  matriz <- matriz[matriz[,1]==0,]
  altura <- 1000*(nrow(matriz)/20)
  fileName <- paste(repoPath,'/results/Figure_4/',organism,'_satStressEnz.png',sep='')
  png(fileName,width=850, height=altura)
  p <- pheatmap(matriz,color = cividis(2),cluster_cols = F,cluster_rows = T, show_rownames = TRUE,scale='none',fontsize = 20)
  dev.off()
}