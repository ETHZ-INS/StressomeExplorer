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

path <- "data/USS-OFT/Output_DLC/"
files <- list.files(path)
model <- load_model_hdf5("data/OFT/ModelOliver.hdf5")
para <- readRDS("data/OFT/ModelOliver_PARAMETERS.Rds")
s2c <- read.table("data/USS-OFT/s2c_OFT.csv",sep = ";", header = T)

pipeline <- function(path){
  Tracking <- ReadDLCDataFromCSV(path, fps = 25)
  Tracking$data$centre <- NULL
  Tracking <- CutTrackingData(Tracking,start = 300,end = 300)
  Tracking <- CalibrateTrackingData(Tracking, "area",in.metric = 42*42, c("tr","tl","bl","br"))
  Tracking <- AddOFTZones(Tracking)
  Tracking <- CleanTrackingData(Tracking, likelihoodcutoff = 0.95, existence.pol = ScalePolygon(Tracking$zones$arena,1.3))
  Tracking <- CalculateAccelerations(Tracking)
  Tracking <- CreateSkeletonData_OFT_v3(Tracking)
  Tracking <- ZscoreNormalizeFeatures(Tracking,omit =names(Tracking$features)[c(1:11,32:35)],type = "mean")
  Tracking$features <- ScaleFeatures(Tracking$features, select = names(Tracking$features)[1:11], factor = 4)
  Tracking$features <- ScaleFeatures(Tracking$features, select = names(Tracking$features)[32:35], factor = 0.1)
  Tracking <- ClassifyBehaviors(Tracking,model,para)
  Tracking <- OFTAnalysis(Tracking, points = "bodycentre" ,movement_cutoff = 5, integration_period = 5)
  return(Tracking)
}

subdat <- RunPipeline(files,path,pipeline)

T45min <- TwoGroupComparisonReport(subdat[s2c[s2c$TimePoint =="45min","DLCFile"]],s2c[s2c$TimePoint =="45min","Group"], FDR.cutoff = 0.05)
T24h <- TwoGroupComparisonReport(subdat[s2c[s2c$TimePoint =="24h","DLCFile"]],s2c[s2c$TimePoint =="24h","Group"], FDR.cutoff = 0.05)

#Behavior videos to check integrity of classifier
behaviors <- c("Supported","Unsupported")

for(j in 1:nrow(s2c)){
  Tracking <- subdat[[s2c[j,"DLCFile"]]]
  for(i in behaviors){
    CreateExampleVideos(Tracking,"classifications",i, n = 2, lag = 0, min.length = 0.5, video = paste("data\\USS\\vid\\",s2c[j,"ID"],".avi",sep = ""), folder = "Examples", name = paste(s2c[j,"ID"],"_B",i,sep = ""))
  }
}

for(i in behaviors){
  vids <- list.files("Examples/")
  get <- paste("_B",i,"_",sep = "")
  write_delim(data.frame("file",paste("'Examples\\",vids[grep(get,vids)],"'",sep = "")), "inputs.txt", col_names = F)
  system(paste("ffmpeg -f concat -safe 0 -i inputs.txt -vcodec copy -acodec copy C",i,".mp4", sep = ""))
}

#Behavior videos to check integrity of classifier
behaviors <- c("Unsupported")

#select 45min swim only
#s2c <- s2c[s2c$TimePoint == "45min" & s2c$Group == "Swim",]

for(j in 1:nrow(s2c)){
  Tracking <- subdat[[s2c[j,"DLCFile"]]]
  for(i in behaviors){
    CreateExampleVideos(Tracking,"classifications",i, n = 40, lag = 0, min.length = 0.25, video = paste("data\\USS-OFT\\vid\\",s2c[j,"ID"],".avi",sep = ""), folder = "Examples", name = paste(s2c[j,"ID"],"_B",i,sep = ""))
  }
}

for(i in behaviors){
  vids <- list.files("Examples/")
  get <- paste("_B",i,"_",sep = "")
  write_delim(data.frame("file",paste("'Examples\\",vids[grep(get,vids)],"'",sep = "")), "inputs.txt", col_names = F)
  system(paste("ffmpeg -f concat -safe 0 -i inputs.txt -vcodec copy -acodec copy C",i,".mp4", sep = ""))
}


