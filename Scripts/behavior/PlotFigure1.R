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
library(stringr)

#barplot OFT

s2c <- read.table("data/USS-OFT/s2c_OFT.csv",sep = ";", header = T)

Results <- read.table("Analyses/USS-OFT/OFT_Unsupervised.csv", header = T, sep = ";")

Resdf <- gather(Results,key = type, value = "value", -"file")
Resdf$is.unsupervsied <- grepl("unsupervised", Resdf$type, fixed = TRUE)
Resdf$is.time <- grepl("time", Resdf$type, fixed = TRUE)

Resdf$origtype <- Resdf$type
Resdf$type <- str_replace(Resdf$type, ".time", "")
Resdf$type <- str_replace(Resdf$type, ".count", "")
Resdf$type <- str_replace(Resdf$type, "classifications.", "")
Resdf$type <- str_replace(Resdf$type, "unsupervised.20.", "Cluster ")
Resdf$type <- str_replace(Resdf$type, "bodycentre.", "")
assign <- match(Resdf$file,s2c$DLCFile)
Resdf$Group <- s2c[assign,"Group"]
Resdf$Animal <- s2c[assign,"Animal"]
Resdf$TimePoint  <- s2c[assign,"TimePoint"]
Resdf$readout <- ifelse(Resdf$is.time,"time","count")

Resdfcount <- Resdf[!Resdf$is.time & Resdf$is.unsupervsied & Resdf$Animal != "USS_12",]
Corder <- paste("Cluster",1:20,sep = " ")
Resdfcount$type <- factor(Resdfcount$type, levels = Corder)
Resdfcount$value <- na.replace(Resdfcount$value)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      sem =  sd(x[[col]], na.rm=TRUE) / sqrt(length(x)))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(Resdfcount, varname="value", 
                    groupnames=c("type", "Group","TimePoint"))

t.tests <- list()
Stats <- NULL
for(i in unique(Resdfcount$type)){
  for(j in unique(Resdfcount$TimePoint)){
    swim <- Resdfcount[Resdfcount$type == i & Resdfcount$TimePoint == j & Resdfcount$Group == "Swim",]
    ctrl <- Resdfcount[Resdfcount$type == i & Resdfcount$TimePoint == j & Resdfcount$Group == "Control",]
    test <- t.test(swim$value,ctrl$value)
    t.tests[[i]] <- test
    Stats <- rbind(Stats, data.frame(type = i, 
                                     TimePoint = j, 
                                     p.value = test$p.value, 
                                     mean.swim = test$estimate[1], 
                                     mean.ctrl = test$estimate[2], 
                                     ratio = test$estimate[1] / test$estimate[2]))
  }
}

p<- ggplot(df2, aes(x=type, y=value, fill=Group, width = 0.75)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values = c("grey","white")) +
  geom_errorbar(aes(ymin=value-sem, ymax=value+sem), width=.2,
                position=position_dodge(.8)) +
  geom_point(data = Resdfcount, aes(type, value, color = Group), position = position_jitterdodge(.3), color = "black") +
  facet_grid(TimePoint~.) + theme_bw()
print(p)
#Export 10x10 inch


#barplot EPM

s2c <- read.table("data/USS-EPM/s2c_EPM.csv",sep = ";", header = T)

Results <- read.table("Analyses/USS-EPM/USS-EPM_Unsupervised_Data.csv", header = T, sep = ";")

Resdf <- gather(Results,key = type, value = "value", -"file")
Resdf$is.unsupervsied <- grepl("unsupervised", Resdf$type, fixed = TRUE)
Resdf$is.time <- grepl("time", Resdf$type, fixed = TRUE)

Resdf$origtype <- Resdf$type
Resdf$type <- str_replace(Resdf$type, ".time", "")
Resdf$type <- str_replace(Resdf$type, ".count", "")
Resdf$type <- str_replace(Resdf$type, "unsupervised.", "Cluster ")
Resdf$type <- str_replace(Resdf$type, "bodycentre.", "")
assign <- match(Resdf$file,s2c$DLCFile)
Resdf$Group <- s2c[assign,"Group"]
Resdf$Animal <- s2c[assign,"Animal"]
Resdf$readout <- ifelse(Resdf$is.time,"time","count")

Resdfcount <- Resdf[Resdf$Animal != "USS_12" & !Resdf$is.time & Resdf$is.unsupervsied & !(Resdf$type %in% c("Cluster 6","Cluster 8","Cluster 9","Cluster 19")),]
Corder <- paste("Cluster",1:20,sep = " ")
Resdfcount$type <- factor(Resdfcount$type, levels = Corder)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      sem =  sd(x[[col]], na.rm=TRUE) / sqrt(length(x)))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

Stats <- NULL
for(i in unique(Resdfcount$type)){
    swim <- Resdfcount[Resdfcount$type == i & Resdfcount$Group == "Swim",]
    ctrl <- Resdfcount[Resdfcount$type == i & Resdfcount$Group == "Control",]
    test <- t.test(swim$value,ctrl$value)
    Stats <- rbind(Stats, data.frame(type = i, 
                                     p.value = test$p.value, 
                                     mean.swim = test$estimate[1], 
                                     mean.ctrl = test$estimate[2], 
                                     ratio = test$estimate[1] / test$estimate[2]))
}

df2 <- data_summary(Resdfcount, varname="value", 
                    groupnames=c("type", "Group"))

p<- ggplot(df2, aes(x=type, y=value, fill=Group, width = 0.75)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values = c("grey","white")) +
  geom_errorbar(aes(ymin=value-sem, ymax=value+sem), width=.2,
                position=position_dodge(.8)) +
  geom_point(data = Resdfcount, aes(type, value, color = Group), position = position_jitterdodge(.3), color = "black")+ theme_bw()
print(p)
#Export 5.5x10 inch


#barplot OFT classical

s2c <- read.table("data/USS-OFT/s2c_OFT.csv",sep = ";", header = T)

Results <- read.table("Analyses/USS-OFT/OFT_Unsupervised.csv", header = T, sep = ";")

Resdf <- gather(Results,key = type, value = "value", -"file")
Resdf$is.unsupervsied <- grepl("unsupervised", Resdf$type, fixed = TRUE)
Resdf$is.time <- grepl("time", Resdf$type, fixed = TRUE)

Resdf$origtype <- Resdf$type
Resdf$type <- str_replace(Resdf$type, ".time", "")
Resdf$type <- str_replace(Resdf$type, ".count", "")
Resdf$type <- str_replace(Resdf$type, "classifications.", "")
Resdf$type <- str_replace(Resdf$type, "unsupervised.20.", "Cluster ")
Resdf$type <- str_replace(Resdf$type, "bodycentre.", "")
assign <- match(Resdf$file,s2c$DLCFile)
Resdf$Group <- s2c[assign,"Group"]
Resdf$Animal <- s2c[assign,"Animal"]
Resdf$TimePoint  <- s2c[assign,"TimePoint"]
Resdf$readout <- ifelse(Resdf$is.time,"time","count")


Resdfcount <- Resdf[Resdf$Animal != "USS_12" & Resdf$origtype %in% c("classifications.Unsupported.count","classifications.Supported.count","bodycentre.distance.moving","bodycentre.center.total.time"),]

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      sem =  sd(x[[col]], na.rm=TRUE) / sqrt(length(x)))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


t.tests <- list()
Stats <- NULL
for(i in unique(Resdfcount$type)){
  for(j in unique(Resdfcount$TimePoint)){
    swim <- Resdfcount[Resdfcount$type == i & Resdfcount$TimePoint == j & Resdfcount$Group == "Swim",]
    ctrl <- Resdfcount[Resdfcount$type == i & Resdfcount$TimePoint == j & Resdfcount$Group == "Control",]
    test <- t.test(swim$value,ctrl$value)
    t.tests[[paste(i,j)]] <- test
    Stats <- rbind(Stats, data.frame(type = i, 
                                     TimePoint = j, 
                                     p.value = test$p.value, 
                                     mean.swim = test$estimate[1], 
                                     mean.ctrl = test$estimate[2], 
                                     ratio = test$estimate[1] / test$estimate[2]))
  }
}
df2 <- data_summary(Resdfcount, varname="value", 
                    groupnames=c("type", "Group","TimePoint"))

p<- ggplot(df2, aes(x=Group, y=value, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values = c("grey","white")) +
  geom_errorbar(aes(ymin=value-sem, ymax=value+sem), width=.2,
                position=position_dodge(.9)) +
  geom_point(data = Resdfcount, aes(Group, value, color = Group), position = position_jitterdodge(), color = "black") +
  facet_grid(type~TimePoint, scales = "free_y") + theme_bw()
print(p)
#Export 10x4 inch


#barplot EPM classical

s2c <- read.table("data/USS-EPM/s2c_EPM.csv",sep = ";", header = T)

Results <- read.table("Analyses/USS-EPM/USS-EPM_Unsupervised_Data.csv", header = T, sep = ";")

Resdf <- gather(Results,key = type, value = "value", -"file")
Resdf$is.unsupervsied <- grepl("unsupervised", Resdf$type, fixed = TRUE)
Resdf$is.time <- grepl("time", Resdf$type, fixed = TRUE)

Resdf$origtype <- Resdf$type
Resdf$type <- str_replace(Resdf$type, ".time", "")
Resdf$type <- str_replace(Resdf$type, ".count", "")
Resdf$type <- str_replace(Resdf$type, "unsupervised.", "Cluster ")
Resdf$type <- str_replace(Resdf$type, "bodycentre.", "")
assign <- match(Resdf$file,s2c$DLCFile)
Resdf$Group <- s2c[assign,"Group"]
Resdf$Animal <- s2c[assign,"Animal"]
Resdf$TimePoint  <- s2c[assign,"TimePoint"]
Resdf$readout <- ifelse(Resdf$is.time,"time","count")


Resdfcount <- Resdf[Resdf$Animal != "USS_12" & Resdf$origtype %in% c("nose.dip","bodycentre.distance.moving","bodycentre.open.total.time","bodycentre.closed.total.time"),]

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      sem =  sd(x[[col]], na.rm=TRUE) / sqrt(length(x)))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

Stats <- NULL
for(i in unique(Resdfcount$type)){
  swim <- Resdfcount[Resdfcount$type == i & Resdfcount$Group == "Swim",]
  ctrl <- Resdfcount[Resdfcount$type == i & Resdfcount$Group == "Control",]
  test <- t.test(swim$value,ctrl$value)
  Stats <- rbind(Stats, data.frame(type = i, 
                                   p.value = test$p.value, 
                                   mean.swim = test$estimate[1], 
                                   mean.ctrl = test$estimate[2], 
                                   ratio = test$estimate[1] / test$estimate[2]))
}

df2 <- data_summary(Resdfcount, varname="value", 
                    groupnames=c("type", "Group"))

p<- ggplot(df2, aes(x=Group, y=value, fill=Group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_fill_manual(values = c("grey","white")) +
  geom_errorbar(aes(ymin=value-sem, ymax=value+sem), width=.2,
                position=position_dodge(.9)) +
  geom_point(data = Resdfcount, aes(Group, value, color = Group), position = position_jitterdodge(), color = "black") +
  facet_grid(type~., scales = "free_y") + theme_bw()

pdf("test.pdf", width = 3, height = 10)
print(p)
dev.off()
#Export 10x4 inch

#TSNE of OFT

subdat <- readRDS("Analyses/USS-OFT/USS_unsupervised.Rds") 
s2c <- read.table("data/USS-OFT/s2c_OFT.csv",sep = ";", header = T)


allx <- list()
ally <- NULL
id <- NULL
frame <- NULL

for(j in names(subdat)){
  allx[[j]] <- CreateTestSet(subdat[[j]], integration_period = 15)$train_x
  ally <- append(ally,na.omit(subdat[[j]]$labels[["unsupervised.20"]]))
  id <- append(id, rep(paste(j),length(subdat[[j]]$frames) -30))
  frame <- append(frame, subdat[[j]]$frames[-c(1:15,(length(subdat[[j]]$frames) - 14):length(subdat[[j]]$frames))])
}

allx <- do.call(rbind, allx)
reducedindex <- sample(1:nrow(allx))[1:(nrow(allx) / 20)]

out <- tsne(t(allx[reducedindex,]),labels = ally[reducedindex],dotsize = 0.5, perplex = 500, seed = 123)
Res2 <- data.frame(out$data, cluster = ally[reducedindex], file = id[reducedindex], frame = frame[reducedindex], index = reducedindex)
assign2 <- match(Res2$file,s2c$DLCFile)
Res2 <- cbind(Res2,s2c[assign2,])

ggplot(Res2, aes(X1,X2, color = cluster)) + 
  geom_point(size = 0.5) + 
  theme_bw()

#Correlation Plot USS-OFT count data

Cordat <- readRDS("Analyses/USS-OFT/USS_unsupervised.Rds")

for(i in names(Cordat)){
  Cordat[[i]]$labels$unsupervised.50 <- NULL
  Cordat[[i]] <- SmoothLabels(Cordat[[i]], 5)
  Cordat[[i]]$Report <- append(Cordat[[i]]$Report,LabelReport(Cordat[[i]])[-c(1:6)])
}

dat <- MultiFileReport(Cordat)
incl <- paste("unsupervised.20.",c(1,2,4:19),".count", sep = "")
incl <- append(incl,c("classifications.Unsupported.count","classifications.Supported.count"))

corrplot(cor(as.matrix(na_replace(dat[,incl]))),title = "Correlation Plot", 
         method = "square", 
         outline = T, 
         addgrid.col = "darkgray", 
         order="hclust")

