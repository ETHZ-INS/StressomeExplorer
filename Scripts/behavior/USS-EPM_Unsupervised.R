setwd("C:/Users/Hackerman/Documents/Lukas/DLCAnalyzer-master/")

source("R/DLCAnalyzer_Functions_v3.R")
library(sp)         #tested with v1.3-2
library(imputeTS)   #tested with v2.7
library(ggplot2)    #tested with v3.1.0
library(ggmap)      #tested with v3.0.0
library(data.table) #tested with v1.12.8
library(cowplot)    #tested with v0.9.4
library(corrplot)   #tested with v0.84
library(keras)      #REQUIRES TENSORFLOW INSTALL. tested with v2.2.5.0
library(readr)
library(tidyr)
library(ggfortify)
library(biganalytics)
library(M3C)

path <- "data/USS-EPM/Output_DLC/"
files <- list.files(path)
source("R/DLCAnalyzer_Functions_v3.R")
s2c <- read.table("data/USS-EPM/s2c_EPM.csv",sep = ";", header = T)

pipeline <- function(path){
  Tracking <- ReadDLCDataFromCSV(path, fps = 25)
  Tracking <- CutTrackingData(Tracking,start = 300,end = 300)
  Tracking <- CalibrateTrackingData(Tracking, method = "distance",in.metric = 65.5, points = c("tl","bl"))
  Tracking <- CleanTrackingData(Tracking, likelihoodcutoff = 0.95)
  zoneinfo <- read.table("example/EPM/EPM_zoneinfo.csv", sep = ";", header = T)
  Tracking <- AddZones(Tracking,zoneinfo)
  Tracking <- EPMAnalysis(Tracking, movement_cutoff = 5,integration_period = 5,points = "bodycentre", nosedips = TRUE)
  Tracking <- CalculateAccelerations(Tracking)
  Tracking <- CreateSkeletonData_EPM(Tracking)
  return(Tracking)
}

subdat <- RunPipeline(files,path,pipeline)

# =========== run kmeans manually===========
allx <- list()
id <- NULL
for(j in names(subdat)){
  subdat[[j]] <- CreateTestSet(subdat[[j]], 15)
  allx[[j]] <- subdat[[j]]$train_x
  id <- append(id, rep(paste(j),nrow(subdat[[j]]$train_x)))
  subdat[[j]]$train_x <- NULL
}

allx <- do.call(rbind, allx)

for(i in 1:ncol(allx)){
  x <- allx[,i]
  allx[,i] <- (x - mean(x))/(sd(x))
}

test <- bigkmeans(allx,centers = 20)

for(j in names(subdat)){
  clust <- test$cluster[id == j]
  subdat[[j]]$labels$unsupervised <- c(rep(NA,subdat[[j]]$ml_integration),as.character(clust),rep(NA,subdat[[j]]$ml_integration))
}

rm(allx)
# =========== end run kmeans manually===========

for(i in names(subdat)){
  subdat[[i]]$Report <- append(subdat[[i]]$Report, LabelReport(subdat[[i]],5))
}
saveRDS(subdat,"USS-EPM_Unsupervised.Rds")

Rep <- TwoGroupComparisonReport(subdat[s2c[,"DLCFile"]],s2c[,"Group"], FDR.cutoff = 0.5)

clusters <- c(1:20)

for(j in 1:nrow(s2c)){
  Tracking <- subdat[[s2c[j,"DLCFile"]]]
  for(i in clusters){
    CreateExampleVideos(Tracking,"unsupervised",i, n = 2, lag = 0, min.length = 0.5, video = paste("data\\USS-EPM\\vid\\",s2c[j,"ID"],".mp4",sep = ""), folder = "Examples", name = paste(s2c[j,"ID"],"_C",i,sep = ""))
  }
}

for(i in clusters){
  vids <- list.files("Examples/")
  get <- paste("_C",i,"_",sep = "")
  write_delim(data.frame("file",paste("'Examples\\",vids[grep(get,vids)],"'",sep = "")), "inputs.txt", col_names = F)
  system(paste("ffmpeg -f concat -safe 0 -i inputs.txt -vcodec copy -acodec copy C",i,".mp4", sep = ""))
}

write.table(MultiFileReport(subdat), "USS-EPM_Unsupervised_Data.csv", sep = ";", row.names = F)
write.table(Rep$Results, "USS-EPM_Unsupervised_Stats.csv", sep = ";", row.names = F)

pdf("USS-EPM_Correlations.pdf", width = 15, height = 15)
print(Rep$PlotsAll)
print(Rep$PlotsSignificant)
CorrelationPlotLabels(subdat, hclust = TRUE)
dev.off()
