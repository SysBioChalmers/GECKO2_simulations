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
repoPath <- getwd()
repoPath <- gsub('/code/analyse_protModels','',repoPath)
setwd(repoPath)
fileName <- paste(repoPath,'/data/SingleCopyOG_All.txt',sep='')
OGDF  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
#convert Kma IDs to provide compatibility with model ids
fileName <- paste(repoPath,'/data/kma_strains_orthologs.txt',sep='')
kmaOGs  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
idxs    <- match(kmaOGs$genes.CBS6556,OGDF$K.marx)
OGDF    <- OGDF[idxs[!is.na(idxs)],]
idxs    <- match(OGDF$K.marx,kmaOGs$genes.CBS6556)
OGDF$K.marx <- kmaOGs$genes.DMKU[idxs]
#convert yli IDs to provide compatibility with model ids
fileName <- paste(repoPath,'/data/yli_strains_orthologs.txt',sep='')
yliOGs  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
idxs    <- match(yliOGs$W29,OGDF$Y.lipo)
OGDF    <- OGDF[idxs[!is.na(idxs)],]
idxs    <- match(OGDF$Y.lipo,yliOGs$W29)
OGDF$Y.lipo <- yliOGs$CLIB122[idxs]
OGDF$Y.lipo <- gsub('YALI0_','YALI0',OGDF$Y.lipo)
fileName <- paste(repoPath,'/data/singleCopyOG.txt',sep='')
write.table(OGDF, file = fileName, row.names = FALSE,quote= FALSE,sep = '\t')

for (organism in c('yli')){
fileName <- paste(repoPath,'/data/',organism,'_enzymeTable.txt',sep='')
enzInfo  <- read.delim(fileName, header = TRUE, sep = "\t",stringsAsFactors=FALSE, na.strings = "NA")
idxs  <- match(OGDF$Y.lipo,enzInfo$genes)
idxs2 <- which(!is.na(idxs))
idxs  <- idxs[!is.na(idxs)]
enzInfo$OG <- rep(NA,nrow(enzInfo))
enzInfo$OG[idxs] <- OGDF$OG..[idxs2]
fileName <- paste(repoPath,'/data/',organism,'_enzymeTable.txt',sep='')
write.table(enzInfo, file = fileName, row.names = FALSE,quote= FALSE,sep = '\t')
}



