library(dplyr)
library(tidyr)
library(ggfortify)
library(ggplot2)
library(varhandle)
library(treemapify)
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
repoPath <- getwd()
repoPath <- gsub('/code/analyse_protModels','',repoPath)
setwd(repoPath)
organism <- 'sce'
fileName <- paste(repoPath,'/results/',organism,'_absUsage.txt',sep='')
dataset  <- read.csv(file = fileName, sep = '\t', header = TRUE, stringsAsFactors = FALSE,dec='.')
#discard zero rows
dataset <- dataset[rowSums(dataset[,5:ncol(dataset)])>0,]
#Load enzymes info
fileName <- paste(repoPath,'/data/',organism,'_enzymeTable.txt',sep='')
enzInfo  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
#Convert absolute enzyme usage data to a mass basis [mmol/gDw] -> [g/gDw]
  proteins <- dataset$enzymes
  avgMW    <- mean(enzInfo$MWs)
  for (i in 1:length(proteins)){
    protein <- proteins[i]
    index   <- grep(protein,enzInfo$enzymes)
    MW      <- avgMW
    if (length(index)>0){
      MW <- enzInfo$MWs[index]
    }
    dataset[i,5:ncol(dataset)] <- dataset[i,5:ncol(dataset)]*MW
  }
  
  #Load KEGG pathways
  fileName  <- paste(repoPath,'/data/keggPathways.txt',sep='')
  KEGG      <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  metGroups <- unique(KEGG$code)
  columnNames <- colnames(dataset)[5:ncol(dataset)]
  proteinAllocation <- data.frame()
  for (metGroup in metGroups){
    dataRows  <- c()
    KEGG_rows <- grep(metGroup,KEGG$code)
    KEGGcodes <- KEGG$KEGGcode[KEGG_rows]
    pathWays  <- KEGG$pathway[KEGG_rows]
    #Add data rows pathway by pathway
    Ptot_path <- c()
    for (i in 1:length(KEGG_rows)){
      #search individual KEGG pathway code in the subSystems array
      kCode    <- KEGGcodes[i]
      pathW    <- pathWays[i]
      print(pathW)
      indexes  <- grep(tolower(pathW),tolower(dataset$subSystems))
      print(length(indexes))
      dataRows <- c(dataRows,indexes)
      if (length(indexes)==0){
        Ptot_path <- rep(0,length(columnNames))
      }else{
        #If KEGG pathway was found in subSystems then get the total protein 
        #burden of the pathway
        Ptot_path <- as.numeric(colSums(dataset[indexes,5:ncol(dataset)]))
      }
      vector <- data.frame(t(Ptot_path))
      colnames(vector) <- columnNames
      proteinAllocation <- rbind(proteinAllocation,cbind(kCode,pathW,metGroup,Ptot_path,columnNames))
    }
  }   
  #Keep pathways with a non-zero total usage
  proteinAllocation <- proteinAllocation[as.numeric(proteinAllocation$Ptot_path)>0,]
  proteinAllocation <- separate(data = proteinAllocation, col = columnNames, into = c("condition", "model"), sep = "_")
  proteinAllocation$Ptot_path  <- unfactor(proteinAllocation$Ptot_path)
  proteinAllocation$condition <- factor(proteinAllocation$condition,levels=c('Std','HiT','LpH','Osm'))
  conditions <- unique(proteinAllocation$condition)
  models <- unique(proteinAllocation$model)
  cumPtotData <- data.frame()
  for (metGroup in metGroups){
    for (model in models){
      for (condition in conditions){
        data <- proteinAllocation[proteinAllocation$condition == condition & proteinAllocation$model == model & proteinAllocation$metGroup == metGroup,]
        Ptot <- sum(data$Ptot_path)
        cumPtotData <- rbind(cumPtotData,cbind(metGroup,model,condition,Ptot))
      }
    }
  }
  #Discard entries for NUC and CofVit
  tempDF <- proteinAllocation
  proteinAllocation     <- proteinAllocation[proteinAllocation$metGroup!='NUC' & proteinAllocation$metGroup!='CofVit',]
  cumPtotData$condition <- factor(cumPtotData$condition,levels=c('Std','HiT','LpH','Osm'))
  df <- cumPtotData[cumPtotData$model=='ecM',]
  colnames(df)[ncol(df)] <- 'ecM'
  df$ecP <-cumPtotData$Ptot[cumPtotData$model=='ecP']
  df$ecP      <- unfactor(df$ecP)
  df$ecM      <- unfactor(df$ecM)
  p <- ggplot(df, aes(x=ecM, y=ecP,color=condition,shape=metGroup)) +
    geom_point(size=6) + 
    scale_x_continuous(limits=c(0,0.03)) + scale_y_continuous(limits = c(0,0.03)) +
    theme_bw(base_size = 30) + scale_color_manual(values=c(rgb(0.45,0.45,0.45),rgb(0.7,0,0.3),rgb(0.8,0.6,0),rgb(0.1,0,0.9))) +
    geom_abline(intercept = 0, slope = 1,linetype="dashed",color = "black", size=0.5) +
    xlab('Protein burden (ecM) [g prot/gDW]') + ylab('Protein burden (ecP) [g prot/gDw]') +
    labs(colour="Condition",shape='') + scale_shape_manual(values = c(15,16,17,18,8))
  plotName <- paste(repoPath,'/results/Figure_4/',organism,'_pairwise_Ptot_pathway.png',sep='')
  png(plotName,width = 650, height = 600)
  plot(p)
  dev.off()
  
  #Plot protein allocation in a treemap for the specific cell phase
  size <- 1000
  for (model in models){
    for (condition in conditions){
      df       <- tempDF[tempDF$model==model & tempDF$condition==condition,]
      df       <- df[df$Ptot_path>0,]
      df$color <- rep('white',nrow(df))
      df$color[df$metGroup=='CEM']    <- 'light blue'
      df$color[df$metGroup=='Lip']    <- 'orange'
      df$color[df$metGroup=='AA']     <- 'light green'
      df$color[df$metGroup=='CofVit'] <- 'red'
      df$color[df$metGroup=='NUC']    <- 'black'
      treeMapPlot <- ggplot(df, aes(area = Ptot_path,
                                fill = color,label=kCode,subgroup=metGroup)) +
    geom_treemap() +
    geom_treemap_subgroup_border(colour = "white") +
    geom_treemap_text(fontface = "italic",colour = "white",place = "centre",grow = F, reflow = T) +
    geom_treemap_subgroup_text(place = "centre",grow = T,alpha = 0.5,colour = "#FAFAFA",min.size = 0)+ 
    scale_fill_identity()
    #plotName <- paste(repoPath,'/results/Figure_4/',organism,'_',model,'_',condition,'_protAllocationTreeMap.png',sep='') 
    #png(plotName,width = size*((sum(df$Ptot_path)^0.5)/(0.1^0.5)), height = size*((sum(df$Ptot_path)^0.5)/(0.1^0.5)))
    #plot(treeMapPlot)
    #dev.off()
    }
  }
 