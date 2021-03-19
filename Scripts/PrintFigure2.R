
install.packages("gridExtra")

library(ggplot2)
library(gridExtra)
library(aggregation)
library(tidyr)

setwd("P:/Lukas/Analysis/20191104_AllPhosphodata/")

# VolcanoPlots
data_volcano <- read.table("VolcanoplotData.csv", sep = ";", header = T)

data_volcano$Condition <- factor(data_volcano$Condition, levels = c("6min","15min","30min","45min"))

FCcutoff <- 0
fdrCutoff <- 0.05
change <- NULL
for(i in 1:dim(data_volcano)[1]){
  if(data_volcano[i,"logFC"] > FCcutoff & data_volcano[i,"FDR"] <= fdrCutoff){
    change<-append(change,"up-regulated")
  }
  else if(data_volcano[i,"logFC"] < -FCcutoff & data_volcano[i,"FDR"] <= fdrCutoff){
    change<-append(change,"down-regulated")
  }
  else{
    change<-append(change,"not significant")
  }
}
data_volcano$change <- change


# ==============VolcanoPlots ==============

library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}



p_volcanos <- ggplot(data_volcano,aes(logFC,PValue,colour=change))+
  geom_point(size = 1, shape = NA)+
  scale_colour_manual(values = c("blue","grey30","red")) + 
  scale_y_continuous(trans=reverselog_trans(10)) +
  scale_x_continuous(breaks = c(-10,0,10), limits = c(-11,11)) +
  theme_bw() +
  facet_grid(Region~Condition)



#Correlation Plots
data_results <- read.table("DEPResults_AllComparisons.csv", sep = ";", header = T)

data_sig <- data_results[data_results$AnySignificant,]
range = 12

data_sig6 <- data_sig
data_sig6$significance <- apply(data_sig6[,c("Exp1_dHC_SW06_vs_Exp1_dHC_CON_FDR","Exp1_vHC_SW06_vs_Exp1_vHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sig6 <- data_sig6[order(-data_sig6$significance),]
p1 <- ggplot(data_sig6,aes(Exp1_dHC_SW06_vs_Exp1_dHC_CON_ratio,Exp1_vHC_SW06_vs_Exp1_vHC_CON_ratio, color = -log10(data_sig6$significance))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("6 min")+
  scale_x_continuous(name = "logFC dHC", limits = c(-10,10), breaks = c(-10,0,10)) +
  scale_y_continuous(name = "logFC vHC", limits = c(-10,10), breaks = c(-10,0,10)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_sig15 <- data_sig
data_sig15$significance <- apply(data_sig15[,c("Exp2_dHC_SW15_vs_Exp2_dHC_CON_FDR","Exp2_vHC_SW15_vs_Exp2_vHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sig15 <- data_sig15[order(-data_sig15$significance),]
p2 <- ggplot(data_sig15,aes(Exp2_dHC_SW15_vs_Exp2_dHC_CON_ratio,Exp2_vHC_SW15_vs_Exp2_vHC_CON_ratio, color = -log10(data_sig15$significance))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("15 min")+
  scale_x_continuous(name = "logFC dHC", limits = c(-10,10), breaks = c(-10,0,10)) +
  scale_y_continuous(name = "logFC vHC", limits = c(-10,10), breaks = c(-10,0,10)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_sig30 <- data_sig
data_sig30$significance <- apply(data_sig30[,c("Exp2_dHC_SW30_vs_Exp2_dHC_CON_FDR","Exp2_vHC_SW30_vs_Exp2_vHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sig30 <- data_sig30[order(-data_sig30$significance),]
p3 <- ggplot(data_sig30,aes(Exp2_dHC_SW30_vs_Exp2_dHC_CON_ratio,Exp2_vHC_SW30_vs_Exp2_vHC_CON_ratio, color = -log10(data_sig30$significance))) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("30 min")+
  scale_x_continuous(name = "logFC dHC", limits = c(-10,10), breaks = c(-10,0,10)) +
  scale_y_continuous(name = "logFC vHC", limits = c(-10,10), breaks = c(-10,0,10)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_sig45 <- data_sig
data_sig45$significance <- apply(data_sig45[,c("Exp2_dHC_SW45_vs_Exp2_dHC_CON_FDR","Exp2_vHC_SW45_vs_Exp2_vHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sig45 <- data_sig45[order(-data_sig45$significance),]
p4 <- ggplot(data_sig45,aes(Exp2_dHC_SW45_vs_Exp2_dHC_CON_ratio,Exp2_vHC_SW45_vs_Exp2_vHC_CON_ratio, color = -log10(data_sig45$significance))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("45 min")+
  scale_x_continuous(name = "logFC dHC", limits = c(-10,10), breaks = c(-10,0,10)) +
  scale_y_continuous(name = "logFC vHC", limits = c(-10,10), breaks = c(-10,0,10)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))
  

data_sigd6v15 <- data_sig
data_sigd6v15$significance <- apply(data_sigd6v15[,c("Exp1_dHC_SW06_vs_Exp1_dHC_CON_FDR","Exp2_dHC_SW15_vs_Exp2_dHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sigd6v15 <- data_sigd6v15[order(-data_sigd6v15$significance),]
pd6v15 <- ggplot(data_sigd6v15,aes(Exp1_dHC_SW06_vs_Exp1_dHC_CON_ratio,Exp2_dHC_SW15_vs_Exp2_dHC_CON_ratio, color = -log10(data_sigd6v15$significance))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("dHC")+
  scale_x_continuous(name = "logFC 6min", limits = c(-10,10), breaks = c(-10,0,10)) +
  scale_y_continuous(name = "logFC 15min", limits = c(-10,10), breaks = c(-10,0,10)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_sigd15v30 <- data_sig
data_sigd15v30$significance <- apply(data_sigd15v30[,c("Exp2_dHC_SW15_vs_Exp2_dHC_CON_FDR","Exp2_dHC_SW30_vs_Exp2_dHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sigd15v30 <- data_sigd15v30[order(-data_sigd15v30$significance),]
pd15v30 <- ggplot(data_sigd15v30,aes(Exp2_dHC_SW15_vs_Exp2_dHC_CON_ratio,Exp2_dHC_SW30_vs_Exp2_dHC_CON_ratio, color = -log10(data_sigd15v30$significance))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("dHC")+
  scale_x_continuous(name = "logFC 15min", limits = c(-10,10), breaks = c(-10,0,10)) +
  scale_y_continuous(name = "logFC 30min", limits = c(-10,10), breaks = c(-10,0,10)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_sigd30v45 <- data_sig
data_sigd30v45$significance <- apply(data_sigd30v45[,c("Exp2_dHC_SW30_vs_Exp2_dHC_CON_FDR","Exp2_dHC_SW45_vs_Exp2_dHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sigd30v45 <- data_sigd30v45[order(-data_sigd30v45$significance),]
pd30v45 <- ggplot(data_sigd30v45,aes(Exp2_dHC_SW30_vs_Exp2_dHC_CON_ratio,Exp2_dHC_SW45_vs_Exp2_dHC_CON_ratio, color = -log10(data_sigd30v45$significance))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("dHC")+
  scale_x_continuous(name = "logFC 30min", limits = c(-10,10), breaks = c(-10,0,10)) +
  scale_y_continuous(name = "logFC 45min", limits = c(-10,10), breaks = c(-10,0,10)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_sigv6v15 <- data_sig
data_sigv6v15$significance <- apply(data_sigv6v15[,c("Exp1_vHC_SW06_vs_Exp1_vHC_CON_FDR","Exp2_vHC_SW15_vs_Exp2_vHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sigv6v15 <- data_sigv6v15[order(-data_sigv6v15$significance),]
pv6v15 <- ggplot(data_sigv6v15,aes(Exp1_vHC_SW06_vs_Exp1_vHC_CON_ratio,Exp2_vHC_SW15_vs_Exp2_vHC_CON_ratio, color = -log10(data_sigv6v15$significance))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("vHC")+
  scale_x_continuous(name = "logFC 6min", limits = c(-10,10), breaks = c(-10,0,10)) +
  scale_y_continuous(name = "logFC 15min", limits = c(-10,10), breaks = c(-10,0,10)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_sigv15v30 <- data_sig
data_sigv15v30$significance <- apply(data_sigv15v30[,c("Exp2_vHC_SW15_vs_Exp2_vHC_CON_FDR","Exp2_vHC_SW30_vs_Exp2_vHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sigv15v30 <- data_sigv15v30[order(-data_sigv15v30$significance),]
pv15v30 <- ggplot(data_sigv15v30,aes(Exp2_vHC_SW15_vs_Exp2_vHC_CON_ratio,Exp2_vHC_SW30_vs_Exp2_vHC_CON_ratio, color = -log10(data_sigv15v30$significance))) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("vHC")+
  scale_x_continuous(name = "logFC 15min", limits = c(-10,10), breaks = c(-10,0,10)) +
  scale_y_continuous(name = "logFC 30min", limits = c(-10,10), breaks = c(-10,0,10)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

data_sigv30v45 <- data_sig
data_sigv30v45$significance <- apply(data_sigv30v45[,c("Exp2_vHC_SW30_vs_Exp2_vHC_CON_FDR","Exp2_vHC_SW45_vs_Exp2_vHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sigv30v45 <- data_sigv30v45[order(-data_sigv30v45$significance),]
pv30v45 <- ggplot(data_sigv30v45,aes(Exp2_vHC_SW30_vs_Exp2_vHC_CON_ratio,Exp2_vHC_SW45_vs_Exp2_vHC_CON_ratio, color = -log10(data_sigv30v45$significance))) + 
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_point(size =1) + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/range, 1),
                         limits = c(0,range)) +
  ggtitle("vHC")+
  scale_x_continuous(name = "logFC 30min", limits = c(-10,10), breaks = c(-10,0,10)) +
  scale_y_continuous(name = "logFC 45min", limits = c(-10,10), breaks = c(-10,0,10)) +
  labs(colour="-log10(FDR)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  guides(color = FALSE)

p_corrs <- grid.arrange(p1, p2,p3,p4, pd6v15,pd15v30,pd30v45,pv6v15,pv15v30,pv30v45, widths = c(1,1,1,1.7), layout_matrix = rbind(c(1, 2, 3, 4),
                                           c(5, 6, 7,NA),c(8, 9, 10,NA)))
#EXPORT AS 8 x 6 inch

grid.arrange(p_volcanos,p_corrs, heights = c(1,0.2,1.3), layout_matrix = rbind(c(1),c(NA),c(2)))
#EXPORT AS 8 x 11 inch

data_profiles <- read.table("Clusterdat_n5.csv", sep = ";", header = T)
data_profiles <- clustering_dat

data_profiles_plot <- NULL
for(i in unique(data_profiles$cluster)){
  dat <- data_profiles[data_profiles$cluster == i,c("log2FC","TimePointNumeric","Region")]
  dat <- as.data.frame(aggregate(log2FC ~ TimePointNumeric + Region, dat, function(x) c(mean = mean(x), N_pep = length(x), sd= sd(x), se = sd(x)/sqrt(length(x)))))
  data_profiles_plot <- rbind(data_profiles_plot,data.frame(dat, cluster = i))
}
data_profiles_plot$log2FC_mean <- data_profiles_plot$log2FC[,"mean"]
data_profiles_plot$log2FC_sd <- data_profiles_plot$log2FC[,"sd"]
data_profiles_plot$log2FC_se <- data_profiles_plot$log2FC[,"se"]
data_profiles_plot$N_pep <- data_profiles_plot$log2FC[,"N_pep"]

data_profiles_plot$SEcolor <- ifelse(data_profiles_plot$Region == "vHC","red","blue")
data_profiles_plot$Title <- apply(data_profiles_plot,1,FUN = function(x){paste("Profile ",x[7],"\n","N = ",x[11], sep = "")})
data_profiles_plot$Title <- as.factor(data_profiles_plot$Title)
data_profiles_plot$Title <- factor(data_profiles_plot$Title, levels = unique(data_profiles_plot[order(-data_profiles_plot$N_pep),"Title"]))

data_profiles_plot <- data_profiles_plot[data_profiles_plot$N_pep > 10,]
pdf("test2.pdf", width = 1.8, height = 5)
p_profiles <- ggplot(data_profiles_plot, aes(TimePointNumeric,log2FC_mean, group = Region)) + 
  geom_hline(yintercept = 0) +
  geom_path(color = data_profiles_plot$SEcolor) +
  scale_x_continuous(breaks = c(6,15,30,45)) +
  xlab( "Time / min") +
  ylab( "mean logFC") +
  geom_ribbon(aes(ymin=data_profiles_plot$log2FC_mean - data_profiles_plot$log2FC_sd, ymax=data_profiles_plot$log2FC_mean + data_profiles_plot$log2FC_sd), color = "grey80", fill = alpha(data_profiles_plot$SEcolor,0.15), size = 0.2) + 
  geom_ribbon(aes(ymin=data_profiles_plot$log2FC_mean - data_profiles_plot$log2FC_se, ymax=data_profiles_plot$log2FC_mean + data_profiles_plot$log2FC_se), color = "black", fill = alpha(data_profiles_plot$SEcolor,1),size = NA) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  facet_grid(Title~., scales = "free")
p_profiles
dev.off()

grid.arrange(p_volcanos,p_corrs, heights = c(1,0.1,1.3), layout_matrix = rbind(c(1),c(NA),c(2)))

ggplot(clustering_dat[clustering_dat$name == "_PVAGGPGAPPAARPPASPSPQR_",], aes(TimePointNumeric,log2FC,color = Region)) + geom_line() + theme_bw()


StatRes <- NULL

comps <- rbind(data.frame(x = 'Exp1_dHC_SW06_vs_Exp1_dHC_CON_ratio',y = 'Exp1_vHC_SW06_vs_Exp1_vHC_CON_ratio'),
               data.frame(x = 'Exp2_dHC_SW15_vs_Exp2_dHC_CON_ratio',y = 'Exp2_vHC_SW15_vs_Exp2_vHC_CON_ratio'),
               data.frame(x = 'Exp2_dHC_SW30_vs_Exp2_dHC_CON_ratio',y = 'Exp2_vHC_SW30_vs_Exp2_vHC_CON_ratio'),
               data.frame(x = 'Exp2_dHC_SW45_vs_Exp2_dHC_CON_ratio',y = 'Exp2_vHC_SW45_vs_Exp2_vHC_CON_ratio'),
               data.frame(x = 'Exp1_vHC_SW06_vs_Exp1_vHC_CON_ratio',y = 'Exp2_vHC_SW15_vs_Exp2_vHC_CON_ratio'),
               data.frame(x = 'Exp2_vHC_SW15_vs_Exp2_vHC_CON_ratio',y = 'Exp2_vHC_SW30_vs_Exp2_vHC_CON_ratio'),
               data.frame(x = 'Exp2_vHC_SW30_vs_Exp2_vHC_CON_ratio',y = 'Exp2_vHC_SW45_vs_Exp2_vHC_CON_ratio'),
               data.frame(x = 'Exp1_dHC_SW06_vs_Exp1_dHC_CON_ratio',y = 'Exp2_dHC_SW15_vs_Exp2_dHC_CON_ratio'),
               data.frame(x = 'Exp2_dHC_SW15_vs_Exp2_dHC_CON_ratio',y = 'Exp2_dHC_SW30_vs_Exp2_dHC_CON_ratio'),
               data.frame(x = 'Exp2_dHC_SW30_vs_Exp2_dHC_CON_ratio',y = 'Exp2_dHC_SW45_vs_Exp2_dHC_CON_ratio'))


for(i in 1:nrow(comps)){
mod <- lm(data_sig[,paste(comps$x[i])]~data_sig[,paste(comps$y[i])])
df <- data.frame(comp = paste(comps$x[i], comps$y[i], sep = " vs. "), slope = mod$coefficients[[2]], intercept = mod$coefficients[[1]], Rsquared = summary(mod)[[8]], pvalue = summary(mod)[[4]][,4][2])
StatRes <- rbind(StatRes,df) 
}
  
write.table(StatRes, "CorrelationStats.csv", sep = ";", row.names = F)
