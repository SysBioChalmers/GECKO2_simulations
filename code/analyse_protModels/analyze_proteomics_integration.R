library(dplyr)
library(fmsb)
organisms  <- c('kma','sce','yli')
orgKeys    <- c('eciSM966','ecYeastGEM','eciYali')
conditions <- c('Std','Hit','LpH','Osm')
if (exists("RStudio.Version")){
  current <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(current)
} else {
  setwd(getSrcDirectory()[1])
}

repoPath <- gsub('/code/analyse_protModels','',current)
methods  <- c('glc')

for (i in 1:length(organisms)){
organism <- organisms[i]
orgKey   <- orgKeys[i]
if (all(organism=='yli')){conditions <- c('Std','Hit','LpH')}

df         <- data.frame(stringsAsFactors = FALSE)
#Load data
fileName <- paste('../',organism,'_scripts/fermentationData.txt',sep='')
expData  <- read.csv(file = fileName, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
dataPath <- paste(repoPath,'/ecModels/',orgKey,'_prot',sep='')
exchData <- data.frame(conditions=conditions,stringsAsFactors=FALSE) 
enzNames <- c()
modProts <- c()

for (j in 1:length(conditions)){
  cond <- conditions[j]
  #Retrieve usages information
  fileName      <- paste(dataPath,'/prot_counts_',cond,'.txt',sep='')
  prot_counts   <- read.csv(file = fileName, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  fileName      <- paste(dataPath,'/enzymeUsages_',cond,'.txt',sep='')
  usages        <- read.csv(file = fileName, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  fileName      <- paste(dataPath,'/modifiedEnzymes_',cond,'.txt',sep='')
  modifications <- read.csv(file = fileName, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
  #Update data frame
  exchData$mass_coverage[which(exchData$conditions==cond)] <- prot_counts$mass_coverage
  exchData$enzNumber[which(exchData$conditions==cond)]     <- prot_counts$matched_prots#/prot_counts$model_prots
  exchData$zeroUsage[which(exchData$conditions==cond)]     <- length(which(usages$usage<=0.01))#/prot_counts$matched_prots
  exchData$fullUsage[which(exchData$conditions==cond)]     <- length(which(usages$usage>=0.99))#/prot_counts$matched_prots
  exchData$avgUsage[which(exchData$conditions==cond)]      <- mean(usages$usage[usages$usage>0])
  exchData$modified [which(exchData$conditions==cond)]     <- sum(modifications$flex_Mass)#/(Ptot[j]*prot_counts$mass_coverage)
  exchData$filtered [which(exchData$conditions==cond)]     <- prot_counts$filtered_prots#/prot_counts$initial_prots
  enzNames <- c(enzNames,usages$prot_IDs)
  modProts <- c(modProts,modifications$prots)
} 
df <- rbind(df,exchData)

filename <- paste('../../results/',organism,'_proteomicsIntegration_abs.txt',sep='')
write.table(df, file = filename, row.names = FALSE, quote = FALSE,sep = '\t')
}
