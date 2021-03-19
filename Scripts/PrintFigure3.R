.libPaths(append(.libPaths(),"C:/Users/lukasv/Documents/RUserLibs"))


library(ggplot2)
library(gridExtra)
library(aggregation)
library(tidyr)


#Plot Figure 3
data_volcano_seq <- read.table("P:/Lukas/Sequencing/20180813_TimeSeriesFST/IndividualTimepoints.csv", sep = ";", header = T)

data_volcano_seq$TimePoint <- factor(data_volcano_seq$TimePoint, levels = c("Swim_45min","Swim_1h30min","Swim_2h","Swim_3h","Swim_4h"))

FCcutoff <- 0
fdrCutoff <- 0.05
change <- NULL
for(i in 1:dim(data_volcano_seq)[1]){
  if(data_volcano_seq[i,"logFC"] > FCcutoff & data_volcano_seq[i,"adj.P.Val"] <= fdrCutoff){
    change<-append(change,"up-regulated")
  }
  else if(data_volcano_seq[i,"logFC"] < -FCcutoff & data_volcano_seq[i,"adj.P.Val"] <= fdrCutoff){
    change<-append(change,"down-regulated")
  }
  else{
    change<-append(change,"not significant")
  }
}
data_volcano_seq$change <- change

p_volcanos_seq <- ggplot(data_volcano_seq,aes(logFC,P.Value,colour=change))+
  geom_point(size = 1, shape = NA)+
  scale_colour_manual(values = c("blue","grey30","red")) + 
  scale_y_continuous(trans=reverselog_trans(10)) +
  scale_x_continuous(breaks = c(-5,0,5), limits = c(-6,6)) +
  theme_bw() +
  facet_grid(Region~TimePoint)


#Correlation Plots
data_seq <- read.table("P:/Lukas/Sequencing/20180813_TimeSeriesFST/IndividualTimepoints.csv",sep = ";", header = T)
data_seq_sig <- unique(data_seq[data_seq$adj.P.Val < 0.05,"EnsembleID"])

range = 121

data_45 <- data_seq[data_seq$TimePoint == "Swim_45min" & data_seq$EnsembleID %in% data_seq_sig,]
data_45 <-  cbind(spread(data_45[,c("EnsembleID","logFC","Region")],value = logFC, key = Region), 
                  Aggregated.FDR = apply(spread(data_45[,c("EnsembleID","adj.P.Val","Region")],
                                                value = adj.P.Val, key = Region)[,c(2,3)],1,FUN = function(x){fisher(x)}))

data_45 <- data_45[order(-data_45$Aggregated.FDR),]
p1 <- ggplot(data_45,aes(dHC,vHC, color = -log10(data_45$Aggregated.FDR))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("45 min")+
  scale_x_continuous(name = "logFC dHC", limits = c(-5,5), breaks = c(-5,0,5)) +
  scale_y_continuous(name = "logFC vHC", limits = c(-5,5), breaks = c(-5,0,5)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_90 <- data_seq[data_seq$TimePoint == "Swim_1h30min" & data_seq$EnsembleID %in% data_seq_sig,]
data_90 <-  cbind(spread(data_90[,c("EnsembleID","logFC","Region")],value = logFC, key = Region), 
                  Aggregated.FDR = apply(spread(data_90[,c("EnsembleID","adj.P.Val","Region")],
                                                value = adj.P.Val, key = Region)[,c(2,3)],1,FUN = function(x){fisher(x)}))
data_90 <- data_90[order(-data_90$Aggregated.FDR),]
p2 <- ggplot(data_90,aes(dHC,vHC, color = -log10(data_90$Aggregated.FDR))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("1h30min")+
  scale_x_continuous(name = "logFC dHC", limits = c(-5,5), breaks = c(-5,0,5)) +
  scale_y_continuous(name = "logFC vHC", limits = c(-5,5), breaks = c(-5,0,5)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_120 <- data_seq[data_seq$TimePoint == "Swim_2h" & data_seq$EnsembleID %in% data_seq_sig,]
data_120 <-  cbind(spread(data_120[,c("EnsembleID","logFC","Region")],value = logFC, key = Region), 
                   Aggregated.FDR = apply(spread(data_120[,c("EnsembleID","adj.P.Val","Region")],
                                                 value = adj.P.Val, key = Region)[,c(2,3)],1,FUN = function(x){fisher(x)}))
data_120 <- data_120[order(-data_120$Aggregated.FDR),]
p3 <- ggplot(data_120,aes(dHC,vHC, color = -log10(data_120$Aggregated.FDR))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("2h")+
  scale_x_continuous(name = "logFC dHC", limits = c(-5,5), breaks = c(-5,0,5)) +
  scale_y_continuous(name = "logFC vHC", limits = c(-5,5), breaks = c(-5,0,5)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_180 <- data_seq[data_seq$TimePoint == "Swim_3h" & data_seq$EnsembleID %in% data_seq_sig,]
data_180 <-  cbind(spread(data_180[,c("EnsembleID","logFC","Region")],value = logFC, key = Region), 
                   Aggregated.FDR = apply(spread(data_180[,c("EnsembleID","adj.P.Val","Region")],
                                                 value = adj.P.Val, key = Region)[,c(2,3)],1,FUN = function(x){fisher(x)}))
data_180 <- data_180[order(-data_180$Aggregated.FDR),]
p4 <- ggplot(data_180,aes(dHC,vHC, color = -log10(data_180$Aggregated.FDR))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("3h")+
  scale_x_continuous(name = "logFC dHC", limits = c(-5,5), breaks = c(-5,0,5)) +
  scale_y_continuous(name = "logFC vHC", limits = c(-5,5), breaks = c(-5,0,5)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_240 <- data_seq[data_seq$TimePoint == "Swim_4h" & data_seq$EnsembleID %in% data_seq_sig,]
data_240 <-  cbind(spread(data_240[,c("EnsembleID","logFC","Region")],value = logFC, key = Region), 
                   Aggregated.FDR = apply(spread(data_240[,c("EnsembleID","adj.P.Val","Region")],
                                                 value = adj.P.Val, key = Region)[,c(2,3)],1,FUN = function(x){fisher(x)}))
data_240 <- data_240[order(-data_240$Aggregated.FDR),]
p5 <- ggplot(data_240,aes(dHC,vHC, color = -log10(data_240$Aggregated.FDR))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("4h")+
  scale_x_continuous(name = "logFC dHC", limits = c(-5,5), breaks = c(-5,0,5)) +
  scale_y_continuous(name = "logFC vHC", limits = c(-5,5), breaks = c(-5,0,5)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))

data_d45v90 <- data_seq[data_seq$TimePoint %in% c("Swim_45min","Swim_1h30min") & data_seq$Region == "dHC" & data_seq$EnsembleID %in% data_seq_sig,]
data_d45v90 <-  cbind(spread(data_d45v90[,c("EnsembleID","logFC","TimePoint")],value = logFC, key = TimePoint), 
                      Aggregated.FDR = apply(spread(data_d45v90[,c("EnsembleID","adj.P.Val","TimePoint")],
                                                    value = adj.P.Val, key = TimePoint)[,c(2,3)],1,FUN = function(x){fisher(x)}))
data_d45v90 <- data_d45v90[order(-data_d45v90$Aggregated.FDR),]
pd45v90 <- ggplot(data_d45v90,aes(Swim_45min,Swim_1h30min, color = -log10(data_d45v90$Aggregated.FDR))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("dHC")+
  scale_x_continuous(name = "logFC 45min", limits = c(-5,5), breaks = c(-5,0,5)) +
  scale_y_continuous(name = "logFC 1h30min", limits = c(-5,5), breaks = c(-5,0,5)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_d90v120 <- data_seq[data_seq$TimePoint %in% c("Swim_1h30min","Swim_2h") & data_seq$Region == "dHC" & data_seq$EnsembleID %in% data_seq_sig,]
data_d90v120 <-  cbind(spread(data_d90v120[,c("EnsembleID","logFC","TimePoint")],value = logFC, key = TimePoint), 
                       Aggregated.FDR = apply(spread(data_d90v120[,c("EnsembleID","adj.P.Val","TimePoint")],
                                                     value = adj.P.Val, key = TimePoint)[,c(2,3)],1,FUN = function(x){fisher(x)}))
data_d90v120 <- data_d90v120[order(-data_d90v120$Aggregated.FDR),]
pd90v120 <- ggplot(data_d90v120,aes(Swim_1h30min,Swim_2h, color = -log10(data_d90v120$Aggregated.FDR))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("dHC")+
  scale_x_continuous(name = "logFC 1h30", limits = c(-5,5), breaks = c(-5,0,5)) +
  scale_y_continuous(name = "logFC 2h", limits = c(-5,5), breaks = c(-5,0,5)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_d120v180 <- data_seq[data_seq$TimePoint %in% c("Swim_2h","Swim_3h") & data_seq$Region == "dHC" & data_seq$EnsembleID %in% data_seq_sig,]
data_d120v180 <-  cbind(spread(data_d120v180[,c("EnsembleID","logFC","TimePoint")],value = logFC, key = TimePoint), 
                        Aggregated.FDR = apply(spread(data_d120v180[,c("EnsembleID","adj.P.Val","TimePoint")],
                                                      value = adj.P.Val, key = TimePoint)[,c(2,3)],1,FUN = function(x){fisher(x)}))
data_d120v180 <- data_d120v180[order(-data_d120v180$Aggregated.FDR),]
pd120v180 <- ggplot(data_d120v180,aes(Swim_2h,Swim_3h, color = -log10(data_d120v180$Aggregated.FDR))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("dHC")+
  scale_x_continuous(name = "logFC 2h", limits = c(-5,5), breaks = c(-5,0,5)) +
  scale_y_continuous(name = "logFC 3h", limits = c(-5,5), breaks = c(-5,0,5)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_d180v240 <- data_seq[data_seq$TimePoint %in% c("Swim_3h","Swim_4h") & data_seq$Region == "dHC" & data_seq$EnsembleID %in% data_seq_sig,]
data_d180v240 <-  cbind(spread(data_d180v240[,c("EnsembleID","logFC","TimePoint")],value = logFC, key = TimePoint), 
                        Aggregated.FDR = apply(spread(data_d180v240[,c("EnsembleID","adj.P.Val","TimePoint")],
                                                      value = adj.P.Val, key = TimePoint)[,c(2,3)],1,FUN = function(x){fisher(x)}))
data_d180v240 <- data_d180v240[order(-data_d180v240$Aggregated.FDR),]
pd180v240 <- ggplot(data_d180v240,aes(Swim_3h,Swim_4h, color = -log10(data_d180v240$Aggregated.FDR))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("dHC")+
  scale_x_continuous(name = "logFC 3h", limits = c(-5,5), breaks = c(-5,0,5)) +
  scale_y_continuous(name = "logFC 4h", limits = c(-5,5), breaks = c(-5,0,5)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)


data_v45v90 <- data_seq[data_seq$TimePoint %in% c("Swim_45min","Swim_1h30min") & data_seq$Region == "vHC" & data_seq$EnsembleID %in% data_seq_sig,]
data_v45v90 <-  cbind(spread(data_v45v90[,c("EnsembleID","logFC","TimePoint")],value = logFC, key = TimePoint), 
                      Aggregated.FDR = apply(spread(data_v45v90[,c("EnsembleID","adj.P.Val","TimePoint")],
                                                    value = adj.P.Val, key = TimePoint)[,c(2,3)],1,FUN = function(x){fisher(x)}))
data_v45v90 <- data_v45v90[order(-data_v45v90$Aggregated.FDR),]
pv45v90 <- ggplot(data_v45v90,aes(Swim_45min,Swim_1h30min, color = -log10(data_v45v90$Aggregated.FDR))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("vHC")+
  scale_x_continuous(name = "logFC 45min", limits = c(-5,5), breaks = c(-5,0,5)) +
  scale_y_continuous(name = "logFC 1h30min", limits = c(-5,5), breaks = c(-5,0,5)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_v90v120 <- data_seq[data_seq$TimePoint %in% c("Swim_1h30min","Swim_2h") & data_seq$Region == "vHC" & data_seq$EnsembleID %in% data_seq_sig,]
data_v90v120 <-  cbind(spread(data_v90v120[,c("EnsembleID","logFC","TimePoint")],value = logFC, key = TimePoint), 
                       Aggregated.FDR = apply(spread(data_v90v120[,c("EnsembleID","adj.P.Val","TimePoint")],
                                                     value = adj.P.Val, key = TimePoint)[,c(2,3)],1,FUN = function(x){fisher(x)}))
data_v90v120 <- data_v90v120[order(-data_v90v120$Aggregated.FDR),]
pv90v120 <- ggplot(data_v90v120,aes(Swim_1h30min,Swim_2h, color = -log10(data_v90v120$Aggregated.FDR))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("vHC")+
  scale_x_continuous(name = "logFC 1h30", limits = c(-5,5), breaks = c(-5,0,5)) +
  scale_y_continuous(name = "logFC 2h", limits = c(-5,5), breaks = c(-5,0,5)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_v120v180 <- data_seq[data_seq$TimePoint %in% c("Swim_2h","Swim_3h") & data_seq$Region == "vHC" & data_seq$EnsembleID %in% data_seq_sig,]
data_v120v180 <-  cbind(spread(data_v120v180[,c("EnsembleID","logFC","TimePoint")],value = logFC, key = TimePoint), 
                        Aggregated.FDR = apply(spread(data_v120v180[,c("EnsembleID","adj.P.Val","TimePoint")],
                                                      value = adj.P.Val, key = TimePoint)[,c(2,3)],1,FUN = function(x){fisher(x)}))
data_v120v180 <- data_v120v180[order(-data_v120v180$Aggregated.FDR),]
pv120v180 <- ggplot(data_v120v180,aes(Swim_2h,Swim_3h, color = -log10(data_v120v180$Aggregated.FDR))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("vHC")+
  scale_x_continuous(name = "logFC 2h", limits = c(-5,5), breaks = c(-5,0,5)) +
  scale_y_continuous(name = "logFC 3h", limits = c(-5,5), breaks = c(-5,0,5)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_v180v240 <- data_seq[data_seq$TimePoint %in% c("Swim_3h","Swim_4h") & data_seq$Region == "vHC" & data_seq$EnsembleID %in% data_seq_sig,]
data_v180v240 <-  cbind(spread(data_v180v240[,c("EnsembleID","logFC","TimePoint")],value = logFC, key = TimePoint), 
                        Aggregated.FDR = apply(spread(data_v180v240[,c("EnsembleID","adj.P.Val","TimePoint")],
                                                      value = adj.P.Val, key = TimePoint)[,c(2,3)],1,FUN = function(x){fisher(x)}))
data_v180v240 <- data_v180v240[order(-data_v180v240$Aggregated.FDR),]
pv180v240 <- ggplot(data_v180v240,aes(Swim_3h,Swim_4h, color = -log10(data_v180v240$Aggregated.FDR))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("vHC")+
  scale_x_continuous(name = "logFC 3h", limits = c(-5,5), breaks = c(-5,0,5)) +
  scale_y_continuous(name = "logFC 4h", limits = c(-5,5), breaks = c(-5,0,5)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)



p_corrs_seq <- grid.arrange(p1, p2,p3,p4, p5,pd45v90,pd90v120,pd120v180,pd180v240,pv45v90,pv90v120,pv120v180,pv180v240, 
                            widths = c(1,1,1,1,1.7), 
                            layout_matrix = rbind(c(1, 2, 3, 4, 5),
                                                  c(6, 7, 8,9,NA),
                                                  c(10, 11, 12,13,NA)))
#EXPORT AS 10 x 6 inch

grid.arrange(p_volcanos_seq,p_corrs_seq, heights = c(1,0.2,1.3), layout_matrix = rbind(c(1),c(NA),c(2)))
#EXPORT AS 10 x 11 inch



# Profiles

StatResults <- read.table("P:/Lukas/Sequencing/20180813_TimeSeriesFST/IndividualTimepoints.csv", sep =";", header = T)
s2c <- read.table(file.path("P:/Lukas/Sequencing/20180813_TimeSeriesFST/metadata", "TimeSeriesFST_s2c.txt"), header = TRUE, stringsAsFactors=FALSE, na.strings = "<NA>")
StatResults$TimePoint <- factor(StatResults$TimePoint, levels = c("Swim_45min","Swim_1h30min","Swim_2h","Swim_3h","Swim_4h"))

data_seq_sig <- unique(data_seq[data_seq$adj.P.Val < 0.05,"EnsembleID"])
expmat_seq <- spread(data_seq[data_seq$EnsembleID %in% data_seq_sig, c("logFC","EnsembleID","Test")], key = Test, value = logFC) 
row.names(expmat_seq) <- expmat_seq[,1]
expmat_seq <- as.matrix(expmat_seq[,-1])
expmat_seq <- expmat_seq[complete.cases(expmat_seq),]


Results <- kmeans(expmat_seq,centers = 25, nstart = 20)
expmat_seq <- data.frame(expmat_seq, cluster = Results$cluster)
expmat_seq$name <- row.names(expmat_seq)

data_profiles <- StatResults[StatResults$EnsembleID %in% expmat_seq$name, c("logFC","EnsembleID","Region","TimePoint")]
assign <- match(data_profiles$EnsembleID,expmat_seq$name)
data_profiles$cluster <- expmat_seq[assign,"cluster"]
assign <- match(data_profiles$EnsembleID,StatResults$EnsembleID)
data_profiles$ShortName <- StatResults[assign,"ShortName"]

data_profiles$TimePointNumeric <-  sapply(data_profiles$TimePoint, FUN = function(x){switch(x,
                                                                                                Swim_45min = as.numeric(45),
                                                                                                Swim_1h30min = as.numeric(90),
                                                                                                Swim_2h = as.numeric(120),
                                                                                                Swim_3h = as.numeric(180),
                                                                                                Swim_4h = as.numeric(240))})
#write.table(data_profiles,"P:/Lukas/Sequencing/20180813_TimeSeriesFST/TimeSeriesCorrelations/Clusters_combined_top25_start20.csv",row.names = F, sep = ";")
data_profiles <- read.table("P:/Lukas/Sequencing/20180813_TimeSeriesFST/TimeSeriesCorrelations/Clusters_combined_top25_start20.csv",header = T, sep = ";")

data_profiles_plot <- NULL
for(i in unique(data_profiles$cluster)){
  dat <- data_profiles[data_profiles$cluster == i,c("logFC","TimePointNumeric","Region")]
  dat <- as.data.frame(aggregate(logFC ~ TimePointNumeric + Region, dat, function(x) c(mean = mean(x), N_pep = length(x), sd= sd(x), se = sd(x)/sqrt(length(x)))))
  data_profiles_plot <- rbind(data_profiles_plot,data.frame(dat, cluster = i))
}
data_profiles_plot$log2FC_mean <- data_profiles_plot$logFC[,"mean"]
data_profiles_plot$log2FC_sd <- data_profiles_plot$logFC[,"sd"]
data_profiles_plot$log2FC_se <- data_profiles_plot$logFC[,"se"]
data_profiles_plot$N_pep <- data_profiles_plot$logFC[,"N_pep"]

data_profiles_plot$SEcolor <- ifelse(data_profiles_plot$Region == "vHC","red","blue")
data_profiles_plot$Title <- apply(data_profiles_plot,1,FUN = function(x){paste("Profile ",x[7],"\n","N = ",x[11], sep = "")})
data_profiles_plot$Title <- as.factor(data_profiles_plot$Title)
data_profiles_plot$Title <- factor(data_profiles_plot$Title, levels = unique(data_profiles_plot[order(-data_profiles_plot$N_pep),"Title"]))

data_profiles_plot <- data_profiles_plot[data_profiles_plot$N_pep > 10,]
pdf("P:/Lukas/Written Papers_Proposals/Stressome 2.0 paper 2019/ProfilesFig3.pdf", height = 12, width = 2)
p_profiles <- ggplot(data_profiles_plot, aes(TimePointNumeric,log2FC_mean, group = Region)) + 
  geom_hline(yintercept = 0) +
  geom_path(color = data_profiles_plot$SEcolor) +
  scale_x_continuous(breaks = c(45,90,120,180,240)) +
  xlab( "Time / min") +
  ylab( "mean logFC") +
  geom_ribbon(aes(ymin=data_profiles_plot$log2FC_mean - data_profiles_plot$log2FC_sd, ymax=data_profiles_plot$log2FC_mean + data_profiles_plot$log2FC_sd), color = "grey80", fill = alpha(data_profiles_plot$SEcolor,0.15), size = 0.2) + 
  geom_ribbon(aes(ymin=data_profiles_plot$log2FC_mean - data_profiles_plot$log2FC_se, ymax=data_profiles_plot$log2FC_mean + data_profiles_plot$log2FC_se), color = "black", fill = alpha(data_profiles_plot$SEcolor,1),size = NA) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  facet_grid(Title~., scales = "free",)
p_profiles
dev.off()

#Correlation stats
StatRes <- NULL
comps <- rbind(data.frame(x = 'vHC_Swim_45min',y = 'dHC_Swim_45min'),
               data.frame(x = 'vHC_Swim_1h30min',y = 'dHC_Swim_1h30min'),
               data.frame(x = 'vHC_Swim_2h',y = 'dHC_Swim_2h'),
               data.frame(x = 'vHC_Swim_3h',y = 'dHC_Swim_3h'),
               data.frame(x = 'vHC_Swim_4h',y = 'dHC_Swim_4h'),
               data.frame(x = 'vHC_Swim_45min',y = 'vHC_Swim_1h30min'),
               data.frame(x = 'vHC_Swim_1h30min',y = 'vHC_Swim_2h'),
               data.frame(x = 'vHC_Swim_2h',y = 'vHC_Swim_3h'),
               data.frame(x = 'vHC_Swim_3h',y = 'vHC_Swim_4h'),
               data.frame(x = 'dHC_Swim_45min',y = 'dHC_Swim_1h30min'),
               data.frame(x = 'dHC_Swim_1h30min',y = 'dHC_Swim_2h'),
               data.frame(x = 'dHC_Swim_2h',y = 'dHC_Swim_3h'),
               data.frame(x = 'dHC_Swim_3h',y = 'dHC_Swim_4h'))

data_sig <- data_seq[data_seq$EnsembleID %in% data_seq_sig,]
data_sig <- spread(data_sig[,c("EnsembleID","logFC","Test")],value = logFC, key = Test)

for(i in 1:nrow(comps)){
  mod <- lm(data_sig[,paste(comps$x[i])]~data_sig[,paste(comps$y[i])])
  df <- data.frame(comp = paste(comps$x[i], comps$y[i], sep = " vs. "), slope = mod$coefficients[[2]], intercept = mod$coefficients[[1]], Rsquared = summary(mod)[[8]], pvalue = summary(mod)[[4]][,4][2])
  StatRes <- rbind(StatRes,df) 
}

write.table(StatRes, "P:/Lukas/Written Papers_Proposals/Stressome 2.0 paper 2019/CorrelationStats_seq.csv", sep = ";", row.names = F)

