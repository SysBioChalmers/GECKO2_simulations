library(reshape)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(matrixStats)
library(pheatmap)
library(viridis)
library(cluster)
library(ggfortify)

#Get repo directory
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
repoPath <- getwd()
repoPath <- gsub('/code/analyse_protModels','',repoPath)
setwd(repoPath)
organisms <- c('sce','kma','yli')
dataTypes <- c('absUsage','protConcs')
fileName <- paste(repoPath,'/data/ecModels_OGs.txt',sep='')
OGs  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
for (dataType in dataTypes){
  matriz <- data.frame(row.names = OGs$OG)
  orgsRow <- c()
  conditions <- c()
for (organism in organisms){
  if (dataType=='absUsage'){fileName  <- paste(repoPath,'/results/',organism,'_absUsage.txt',sep='')}
  if (dataType=='protConcs'){fileName <- paste(repoPath,'/results/ecModels_',organism,'_protConcs.txt',sep='')}
  protConcs  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
  if (dataType=='absUsage'){protConcs<- protConcs[,-grep('_ecM',colnames(protConcs))] }
  column <- grep(organism,colnames(OGs))
  idxs   <- match(OGs[,column],protConcs$genes)
  numCol <- 5
  if (dataType=='protConcs'){numCol <- 6}
  
  if (all(organism=='sce')){
    df <- data.frame(protConcs[idxs,1:4])
  }
  matriz     <- cbind(matriz,protConcs[idxs,numCol:ncol(protConcs)])
  orgConds   <- colnames(protConcs[numCol:ncol(protConcs)])
  conditions <- c(conditions,orgConds)
  orgsRow    <- c(orgsRow,rep(organism,length(orgConds)))
}
colnames(matriz) <- gsub('_ecP','',colnames(matriz))
conditions <- gsub('_ecP','',conditions)

df$OGs    <- OGs$OG
organisms <- orgsRow
matriz    <- as.data.frame(t(matriz))
matriz <- do.call(data.frame,lapply(matriz, function(x) replace(x, is.infinite(x),0)))
matriz <- do.call(data.frame,lapply(matriz, function(x) replace(x, is.na(x),0)))
matriz$conditions <- factor(conditions,levels=c('Std','HiT','LpH','Osm'))
matriz$organisms  <- factor(organisms,levels=c('sce','kma','yli'))

PCAdata  <- prcomp(matriz[,(1:(ncol(matriz)-2))], center = TRUE,scale = FALSE,retx=TRUE)
p        <- autoplot(PCAdata,data = matriz,shape= 'organisms',colour = 'conditions',size = 6,frame = TRUE, frame.type = 'norm')
p        <- p + theme_bw(base_size = 24)
p        <- p + scale_color_manual(values=c(rgb(0.45,0.45,0.45),rgb(0.7,0,0.3),rgb(0.8,0.5,0),rgb(0.1,0,0.9)))
p        <- p + stat_ellipse(aes(group=organisms),level=0.95)
plotName <- paste(repoPath,'/results/Figure_S3/', dataType ,'_OGs_PCA.png',sep='') 
png(plotName,width = 650, height = 600)
plot(p)
dev.off()
}
