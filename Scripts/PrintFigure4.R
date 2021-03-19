.libPaths(append(.libPaths(),"C:/Users/lukasv/Documents/RUserLibs"))

library(ggplot2)
library(gridExtra)
library(aggregation)
library(tidyr)

data_spl_seq <- read.table("P:/Lukas/Written Papers_Proposals/Stressome 2.0 paper 2019/SplicedvsUnspliced.csv", sep =";", header = T)
clusters <- read.table("P:/Lukas/Written Papers_Proposals/Stressome 2.0 paper 2019/Clusters_combined_top25_start20.csv", sep = ";", header = T)

data_spl_seq_long <- data_spl_seq[,c("logFC.TimePoint45min","logFC.TimePoint1h30min","logFC.TimePoint2h","logFC.TimePoint3h","logFC.TimePoint4h","logFC.typeunspliced.TimePoint45min","logFC.typeunspliced.TimePoint1h30min","logFC.typeunspliced.TimePoint2h","logFC.typeunspliced.TimePoint3h","logFC.typeunspliced.TimePoint4h")]
data_spl_seq_long$ShortName <- row.names(data_spl_seq_long)
data_spl_seq_long <- gather(data_spl_seq_long, key = "Test",value = "logFC", -ShortName)

data_spl_seq_long$TimePointNumeric <-  sapply(data_spl_seq_long$Test, FUN = function(x){switch(x,
                                                                                         logFC.TimePoint45min = as.numeric(45),
                                                                                         logFC.TimePoint1h30min = as.numeric(90),
                                                                                         logFC.TimePoint2h = as.numeric(120),
                                                                                         logFC.TimePoint3h = as.numeric(180),
                                                                                         logFC.TimePoint4h = as.numeric(240),
                                                                                         logFC.typeunspliced.TimePoint45min = as.numeric(45),
                                                                                         logFC.typeunspliced.TimePoint1h30min = as.numeric(90),
                                                                                         logFC.typeunspliced.TimePoint2h = as.numeric(120),
                                                                                         logFC.typeunspliced.TimePoint3h = as.numeric(180),
                                                                                         logFC.typeunspliced.TimePoint4h = as.numeric(240))})
data_spl_seq_long$type <-  sapply(data_spl_seq_long$Test, FUN = function(x){switch(x,
                                                                                 logFC.TimePoint45min = "spliced",
                                                                                 logFC.TimePoint1h30min = "spliced",
                                                                                 logFC.TimePoint2h = "spliced",
                                                                                 logFC.TimePoint3h = "spliced",
                                                                                 logFC.TimePoint4h = "spliced",
                                                                                 logFC.typeunspliced.TimePoint45min = "unspliced",
                                                                                 logFC.typeunspliced.TimePoint1h30min = "unspliced",
                                                                                 logFC.typeunspliced.TimePoint2h = "unspliced",
                                                                                 logFC.typeunspliced.TimePoint3h = "unspliced",
                                                                                 logFC.typeunspliced.TimePoint4h = "unspliced")})
unique(clusters[clusters$cluster == 1,"ShortName"])

data_profiles_plot <- NULL
for(i in unique(clusters$cluster)){
  tryCatch({
  dat <- data_spl_seq_long[data_spl_seq_long$ShortName %in% unique(clusters[clusters$cluster == i,"ShortName"]),c("logFC","TimePointNumeric","type")]
  dat <- as.data.frame(aggregate(logFC ~ TimePointNumeric + type, dat, function(x) c(mean = mean(x), N_pep = length(x), sd= sd(x), se = sd(x)/sqrt(length(x)))))
  data_profiles_plot <- rbind(data_profiles_plot,data.frame(dat, cluster = i))
  }, error = function(e){print(e)})
}
data_profiles_plot$log2FC_mean <- data_profiles_plot$logFC[,"mean"]
data_profiles_plot$log2FC_sd <- data_profiles_plot$logFC[,"sd"]
data_profiles_plot$log2FC_se <- data_profiles_plot$logFC[,"se"]
data_profiles_plot$N_pep <- data_profiles_plot$logFC[,"N_pep"]

data_profiles_plot$SEcolor <- ifelse(data_profiles_plot$type == "spliced","red","blue")
data_profiles_plot$Title <- apply(data_profiles_plot,1,FUN = function(x){paste("Profile ",x[7],"\n","N = ",x[11], sep = "")})
data_profiles_plot$Title <- as.factor(data_profiles_plot$Title)
data_profiles_plot$Title <- factor(data_profiles_plot$Title, levels = unique(data_profiles_plot[order(-data_profiles_plot$N_pep),"Title"]))

data_profiles_plot <- data_profiles_plot[data_profiles_plot$N_pep > 10,]
pdf("P:/Lukas/Written Papers_Proposals/Stressome 2.0 paper 2019/ProfilesFig4.pdf", height = 12, width = 2)
p_profiles <- ggplot(data_profiles_plot, aes(TimePointNumeric,log2FC_mean, group = type)) + 
  geom_hline(yintercept = 0) +
  geom_path(color = data_profiles_plot$SEcolor) +
  scale_x_continuous(breaks = c(45,90,120,180,240)) +
  xlab( "Time / min") +
  ylab( "mean logFC") +
  geom_ribbon(aes(ymin=data_profiles_plot$log2FC_mean - data_profiles_plot$log2FC_se, ymax=data_profiles_plot$log2FC_mean + data_profiles_plot$log2FC_se), color = "black", fill = alpha(data_profiles_plot$SEcolor,1),size = NA) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  facet_grid(Title~., scales = "free",)
p_profiles
dev.off()


