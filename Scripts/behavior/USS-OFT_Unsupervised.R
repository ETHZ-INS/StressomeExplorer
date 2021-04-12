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
source("R/DLCAnalyzer_Functions_v3.R")
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
  subdat[[i]]$train_x <- NULL
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
  subdat[[j]]$labels$unsupervised.20 <- c(rep(NA,subdat[[j]]$ml_integration),as.character(clust),rep(NA,subdat[[j]]$ml_integration))
}

test <- bigkmeans(allx,centers = 50)
for(j in names(subdat)){
  clust <- test$cluster[id == j]
  subdat[[j]]$labels$unsupervised.50 <- c(rep(NA,subdat[[j]]$ml_integration),as.character(clust),rep(NA,subdat[[j]]$ml_integration))
}

rm(allx)

saveRDS(subdat,"USS_unsupervised.Rds")
# =========== end run kmeans manually===========


T45min <- subdat[s2c[s2c$TimePoint == "45min","DLCFile"]]
for(i in names(T45min)){
  T45min[[i]]$labels$unsupervised.50 <- NULL
  T45min[[i]] <- OFTAnalysis(T45min[[i]], points = "bodycentre" ,movement_cutoff = 5, integration_period = 5)
}
Rep_45min.20 <- TwoGroupComparisonReport(T45min,by = s2c[s2c$TimePoint == "45min","Group"], FDR.cutoff = 0.5)

T45min <- subdat[s2c[s2c$TimePoint == "45min","DLCFile"]]
for(i in names(T45min)){
  T45min[[i]]$labels$unsupervised.20 <- NULL
  T45min[[i]] <- OFTAnalysis(T45min[[i]], points = "bodycentre" ,movement_cutoff = 5, integration_period = 5)
}
Rep_45min.50 <- TwoGroupComparisonReport(T45min,by = s2c[s2c$TimePoint == "45min","Group"], FDR.cutoff = 0.5)

T45min <- subdat[s2c[s2c$TimePoint == "24h","DLCFile"]]
for(i in names(T45min)){
  T45min[[i]]$labels$unsupervised.50 <- NULL
  T45min[[i]] <- OFTAnalysis(T45min[[i]], points = "bodycentre" ,movement_cutoff = 5, integration_period = 5)
}
Rep_24h.20 <- TwoGroupComparisonReport(T45min,by = s2c[s2c$TimePoint == "24h","Group"], FDR.cutoff = 0.5)

T45min <- subdat[s2c[s2c$TimePoint == "24h","DLCFile"]]
for(i in names(T45min)){
  T45min[[i]]$labels$unsupervised.20 <- NULL
  T45min[[i]] <- OFTAnalysis(T45min[[i]], points = "bodycentre" ,movement_cutoff = 5, integration_period = 5)
}
Rep_24h.50 <- TwoGroupComparisonReport(T45min,by = s2c[s2c$TimePoint == "24h","Group"], FDR.cutoff = 0.5)


clusters <- c(1:20)

for(j in 1:nrow(s2c)){
  Tracking <- subdat[[s2c[j,"DLCFile"]]]
  for(i in clusters){
    CreateExampleVideos(Tracking,"unsupervised.20",i, n = 2, lag = 0, min.length = 0.5, video = paste("data\\USS\\vid\\",s2c[j,"ID"],".avi",sep = ""), folder = "Examples", name = paste(s2c[j,"ID"],"_C",i,sep = ""))
  }
}

for(i in clusters){
  vids <- list.files("Examples/")
  get <- paste("_C",i,"_",sep = "")
  write_delim(data.frame("file",paste("'Examples\\",vids[grep(get,vids)],"'",sep = "")), "inputs.txt", col_names = F)
  system(paste("ffmpeg -f concat -safe 0 -i inputs.txt -vcodec copy -acodec copy C",i,".mp4", sep = ""))
}

Cordat <- subdat
for(i in names(Cordat)){
  Cordat[[i]]$labels$unsupervised.50 <- NULL
  Cordat[[i]] <- SmoothLabels(Cordat[[i]], 5)
  Cordat[[i]]$Report <- append(Cordat[[i]]$Report,LabelReport(Cordat[[i]])[-c(1:6)])
}
  
  
pdf("USS-OFT.pdf", height = 15, width = 15)
print(Rep_24h.20$PlotsAll)
print(Rep_45min.20$PlotsAll)
print(Rep_45min.20$PlotsSignificant)
CorrelationPlotLabels(Cordat,hclust = TRUE)
dev.off()

write.table(MultiFileReport(Cordat), "Analyses/USS-OFT/OFT_Unsupervised.csv", sep = ";", row.names = F)


