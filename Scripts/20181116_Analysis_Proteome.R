.libPaths(append(.libPaths(),"C:/Users/lukasv/Documents/%HOMESHARE%/R3UserLibs"))
.libPaths(append(.libPaths(),"C:/Users/lukasv/Documents/RUserLibs"))


#install.packages("gplots", lib = "C:/Users/lukasv/Documents/%HOMESHARE%/R3UserLibs")
#install.packages("lme4",lib = "C:/Users/lukasv/Documents/%HOMESHARE%/R3UserLibs")
#install.packages("lmerTest",lib = "C:/Users/lukasv/Documents/%HOMESHARE%/R3UserLibs")
library("lme4")
library("lmerTest")
library(tidyr)
library(ggplot2)
library(gplots)
library(edgeR)

#===========================================================================================================================================
#=======================================================LOAD DATA AND ASSING GROUPING VARIABLES=============================================
#===========================================================================================================================================

setwd("P:/Lukas/Experiments/20180221_LvZ1_FST_4H_Proteomics/20181116_Reanalysis")

Proteins <- read.table("Data/20181113_TimeseriesFST_Reanalysis_ProteinReport_expanded.csv", header = T, sep = ",", na.strings = "NaN")

s2c <- read.delim("metadata/s2c.txt",header = T)
blocks <- read.delim("metadata/blocks.txt",header = T)


assign <- match(Proteins$R.Condition,s2c$ConNum)
Proteins <- cbind(Proteins, Condition = s2c[assign,"Condition"], Region = s2c[assign,"Region"], SubRegion = s2c[assign,"SubRegion"])
assign <- match(Proteins$R.Replicate,blocks$Replicate)
Proteins <- cbind(Proteins, Block = blocks[assign,"Block"])
AnimalNumber <- do.call('rbind', strsplit(as.character(Proteins$R.FileName),'_',fixed=TRUE))[,6]
AnimalNumber <- do.call('rbind', strsplit(as.character(AnimalNumber),'v',fixed=TRUE))[,1]
AnimalNumber <- paste("id_",do.call('rbind', strsplit(as.character(AnimalNumber),'d',fixed=TRUE))[,1],sep = "")
Proteins <- cbind(Proteins,AnimalNumber)

Proteins$RawQuant <- Proteins$PG.Quantity
Proteins$PG.Quantity <- log2(Proteins$PG.Quantity)

assign <- match(Peptides$R.Condition,s2c$ConNum)
Peptides <- cbind(Peptides, Condition = s2c[assign,"Condition"], Region = s2c[assign,"Region"], SubRegion = s2c[assign,"SubRegion"])
assign <- match(Peptides$R.Replicate,blocks$Replicate)
Peptides <- cbind(Peptides, Block = blocks[assign,"Block"])

CompleteMappingTable <- Proteins[,c("R.FileName","Condition","Region","SubRegion","Block")] 
CompleteMappingTable <- CompleteMappingTable[!duplicated(CompleteMappingTable),]
AccessionNameMapping <- Proteins[,c("PG.ProteinAccessions","PG.Genes")]
AccessionNameMapping <- AccessionNameMapping[!duplicated(CompleteMappingTable),]




#===========================================================================================================================================
#=======================================================VISUALIZE DATA======================================================================
#===========================================================================================================================================

# ======== Protein Level =============

# Main effect (Swim vs Homecage) in defined Subregion and Region
Region <- "vHC"
SubRegion <- "DG"
# Raw Data
ggplot(Proteins[Proteins$PG.Genes =="Gapdh"& Proteins$SubRegion == SubRegion & Proteins$Region == Region,],aes(PG.Genes,PG.Quantity, color=Condition)) + geom_boxplot(outlier.shape = NA) +  geom_point(position = position_jitterdodge()) 
# With Blocking
ggplot(Proteins[Proteins$PG.Genes =="Rps5"& Proteins$SubRegion == SubRegion & Proteins$Region == Region,],aes(PG.Genes,PG.Quantity, color=Condition)) + geom_boxplot(outlier.shape = NA) + facet_grid(.~Block) +  geom_point(position = position_jitterdodge()) 


# Main effect (Swim vs Homecage) in all RegionsxSubregions
ggplot(Proteins[Proteins$PG.Genes =="Psmd13",],aes(PG.Genes,PG.Quantity, color=Condition)) + geom_boxplot(outlier.shape = NA) + facet_grid(SubRegion~Region) +  geom_point(position = position_jitterdodge()) 
# With Blocking
ggplot(Proteins[Proteins$PG.Genes =="Psmd13",],aes(PG.Genes,PG.Quantity, color=Condition)) + geom_boxplot(outlier.shape = NA) + facet_grid(SubRegion~Region~Block) +  geom_point(position = position_jitterdodge()) 


# vHC vs dHC in Subregions
ggplot(Proteins[Proteins$PG.Genes =="Shisa7",],aes(PG.Genes,PG.Quantity, color=Region)) + geom_boxplot(outlier.shape = NA) + facet_grid(.~SubRegion) +  geom_point(position = position_jitterdodge()) 
# Subregions in Regions
ggplot(Proteins[Proteins$PG.Genes =="Dpp10",],aes(PG.Genes,PG.Quantity, color=SubRegion)) + geom_boxplot(outlier.shape = NA) + facet_grid(.~Region) +  geom_point(position = position_jitterdodge()) 


# ======== Peptide Level =============
ggplot(Peptides[Peptides$PG.Genes =="Sod1" & Peptides$EG.IsProteotypic == "True" & Peptides$SubRegion == SubRegion & Peptides$Region == Region,],aes(EG.PrecursorId,FG.Quantity, color=Condition)) + 
  geom_boxplot(outlier.shape = NA)  + 
  facet_grid(.~Block) + 
  geom_point(position = position_jitterdodge())+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_log10()


#===========================================================================================================================================
#=======================================================HEATMAPS============================================================================
#===========================================================================================================================================
Proteins = cbind(Proteins,SampleName = paste(substring(Proteins$Condition,1,1), "-", substring(Proteins$Region,1,1),Proteins$SubRegion,"-",Proteins$R.Replicate,sep = ""))
ProteinMatrix <- spread(Proteins[,c("PG.ProteinAccessions","SampleName","PG.Quantity")], SampleName, PG.Quantity)
row.names(ProteinMatrix) <- ProteinMatrix[,1]
ProteinMatrix <- as.matrix(ProteinMatrix[,2:dim(ProteinMatrix)[2]])
keep <- rowSums(is.na(ProteinMatrix)) == 0
ProteinMatrix <- ProteinMatrix[keep,]
ProteinMatrix <- t(apply(ProteinMatrix, 1, function(x)((x - mean(x)) / sd(x))))
heatmap.2(ProteinMatrix, 
          trace = "none",col=redgreen(16), 
          breaks = seq(from = -2, to = 2, by = 0.25), 
          key.xlab = "#sd distance from mean",
          xlab = "Sample", ylab = "Protein", labRow = NA)


#===========================================================================================================================================
#=======================================================PRINCIPAL COMPONENT ANALYSIS (PCA)==================================================
#===========================================================================================================================================
ProteinPCA <- prcomp(ProteinMatrix,
                 center = TRUE,
                 scale. = TRUE)
summary(ProteinPCA)

circle <- function(center = c(0, 0), npoints = 100) {
  r = 1
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0, 0), npoints = 100)

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(ProteinMatrix, ProteinPCA$x))

groups <- c("H-dCA1","S-dCA1","H-dCA3","S-dCA3","H-dDG-","S-dDG-","H-vCA1","S-vCA1","H-vCA3","S-vCA3","H-vDG-","S-vDG-")
colors <- c("a","b","c","d","e","f","g","h","i","j","k","l")
customcols <- c("skyblue4","skyblue2","red4","red2","darkolivegreen4","darkolivegreen3","turquoise4","turquoise2","orange3","orange","green3","green1")
colormapping <- data.frame(groups,colors)
assign <- match(substring(rownames(correlations),1,6),colormapping$groups)
textcolor <- colormapping[assign,2]
shape <- textcolor %in% c("a","b","c","d","e","f")
color <- textcolor
sapply(color,FUN = function(x){switch(x){
  
}})
color[color %in% c("a","b","g","h")] = "a"
color[color %in% c("c","d","i","j")] = "c"
color[color %in% c("e","f","k","l")] = "e"
customcols <- c("red","blue","green")


# geom_path will do open circles
ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "black") + 
  geom_text(data = correlations, aes(x = PC1, y = PC3, label = rownames(correlations), color = textcolor)) + scale_color_manual(values = customcols) + theme_bw()+
  geom_hline(yintercept = 0, colour = "black") + geom_vline(xintercept = 0, 
                                                             colour = "black") + xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + labs(x = "PC1 19.9% of variance", 
                                                                                                                           y = "PC3 12.8% of variance") + ggtitle("Circle of correlations") + theme(legend.position="none")

ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "black") + 
  geom_point(data = correlations, aes(x = PC1, y = PC3, shape = shape, color = color)) + scale_color_manual(values = customcols) + theme_bw()+
  geom_hline(yintercept = 0, colour = "black") + geom_vline(xintercept = 0, 
                                                            colour = "black") + xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + labs(x = "PC1 19.9% of variance", 
                                                                                                                         y = "PC3 12.8% of variance") + ggtitle("Circle of correlations") + theme(legend.position="none")

#===========================================================================================================================================
#=======================================================TWO GROUP COMPARISONS=====================================================
#===========================================================================================================================================


# ============== Swim vs Homecage within selected region and sub-region==============:
AllResults <-NULL
for(i in unique(Proteins$SubRegion)){
  for(j in unique(Proteins$Region)){
    
    print(paste(i, j , sep = " "))

SubSet <- Proteins[Proteins$SubRegion == i & Proteins$Region == j,]
group1 = "Swim"
group2 = "Homecage"

# Create Matrix to filter out low abundante Proteins, Create Target list for Testing
FilterMatrix <- spread(SubSet[,c("PG.ProteinAccessions","R.FileName","PG.Quantity")], R.FileName, PG.Quantity)
keep <- rowSums(is.na(FilterMatrix)) <= 4
ProteinAccesions <- FilterMatrix[keep,1]

# Performa ANOVA for all proteins
Results <- NULL
for(k in 1:length(ProteinAccesions)){
  TestSet <- SubSet[SubSet$PG.ProteinAccessions == ProteinAccesions[k],]
  ANOVA <- aov(formula = PG.Quantity ~ Condition + Error(Block), TestSet)
  pvalue <- as.numeric(unlist(summary(ANOVA))["Error: Within.Pr(>F)1"])
  log2FC <- as.numeric(mean(TestSet[TestSet$Condition == group1,]$PG.Quantity, na.rm = T) - mean(TestSet[TestSet$Condition == group2,]$PG.Quantity, na.rm = T))
  Results <- rbind(Results,data.frame(paste(ProteinAccesions[k]),pvalue, log2FC))
}

# Create Results table and perform FDR testing
names(Results) <- c("Accession","Pvalue","Log2FC")
assign <- match(Results$Accession,AccessionNameMapping$PG.ProteinAccessions)
Results <- cbind(Results, ShortName = AccessionNameMapping[assign,"PG.Genes"])
Results <- cbind(Results,fdr = p.adjust(Results$Pvalue,method = "fdr"))
Results <- cbind(Results,Region = j)
Results <- cbind(Results,SubRegion = i)

#Append to Final Result Table
AllResults <- rbind(AllResults, Results)
  }
}

#write final results Table
write.table(AllResults,"StatisticalResults_SwimvsHomecage_AllRegionsSubRegions.csv",sep=";",row.names =F)

#===========================================================================================================================================
#=======================================================VOLCANO PLOTS OF TWO GROUP COMPARISONS==============================================
#===========================================================================================================================================

#==============Determine DE based on FC and FDR cutoff ==============
FCcutoff <- 0
fdrCutoff <- 0.5
AllResults$change <- NULL
change <- NULL
for(i in 1:dim(AllResults)[1]){
  if(AllResults[i,"Log2FC"] > FCcutoff & AllResults[i,"fdr"] <= fdrCutoff){
    change<-append(change,"up-regulated")
  }
  else if(AllResults[i,"Log2FC"] < -FCcutoff & AllResults[i,"fdr"] <= fdrCutoff){
    change<-append(change,"down-regulated")
  }
  else{
    change<-append(change,"not significant")
  }
}
AllResults<-cbind(AllResults,change)


# ==============VolcanoPlots ==============

library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

ggplot(AllResults,aes(Log2FC,Pvalue,colour=change))+
  geom_point()+
  scale_colour_manual(values = c("blue","black","red")) + 
  scale_y_continuous(trans=reverselog_trans(10)) + 
  scale_x_continuous(limits = c(-2,2))+ 
  facet_grid(Region~SubRegion)+
  theme_bw()


ggplot(AllResults,aes(Log2FC,fdr,colour=change))+geom_point()+scale_colour_manual(values = c("blue","black","red"))  + scale_x_continuous(limits = c(-2,2))+ scale_y_reverse(limits = c(1,0)) + facet_grid(Region~SubRegion)



#===========================================================================================================================================
#=======================================================LMER ANALYSIS==============================================
#===========================================================================================================================================


data <- Proteins[Proteins$PG.Genes =="Rps5" & Proteins$SubRegion == "DG" & Proteins$Region == "vHC",]

data$Group <- as.factor(paste(data$Condition,data$Region,data$SubRegion,sep = "_"))
levels(data$Group)

design <- lmer(PG.Quantity ~  Condition + (1|Block) , data = data)


#===========================================================================================================================================
#=======================================================Region / SubRegion Analysis ==============================================
#===========================================================================================================================================

# ==============CA1 vs CA3 / CA1 vs DG / CA3 vs DG for selected regions:============================
AllResults <-NULL
Region <- "dHC"
group1 <- "CA3"
group2 <- "DG"
SubSet <- Proteins[Proteins$Region == Region & (Proteins$SubRegion == group1 | Proteins$SubRegion == group2),]

# Create Matrix to filter out low abundante Proteins, Create Target list for Testing
FilterMatrix <- spread(SubSet[,c("PG.ProteinAccessions","R.FileName","PG.Quantity")], R.FileName, PG.Quantity)
keep <- rowSums(is.na(FilterMatrix)) <= 4
ProteinAccesions <- FilterMatrix[keep,1]

# Performa ANOVA for all proteins
Results <- NULL
for(i in 1:length(ProteinAccesions)){
  TestSet <- SubSet[SubSet$PG.ProteinAccessions == ProteinAccesions[i],]
  ANOVA <- aov(formula = PG.Quantity ~ SubRegion + Error(AnimalNumber), TestSet)
  pvalue <- as.numeric(unlist(summary(ANOVA))["Error: Within.Pr(>F)1"])
  log2FC <- as.numeric(mean(TestSet[TestSet$SubRegion == group1,]$PG.Quantity, na.rm = T) - mean(TestSet[TestSet$SubRegion == group2,]$PG.Quantity, na.rm = T))
  Results <- rbind(Results,data.frame(paste(ProteinAccesions[i]),pvalue, log2FC))
}

# Create Results table and perform FDR testing
names(Results) <- c("Accession","Pvalue","Log2FC")
assign <- match(Results$Accession,AccessionNameMapping$PG.ProteinAccessions)
Results <- cbind(Results, ShortName = AccessionNameMapping[assign,"PG.Genes"])
Results <- cbind(Results,fdr = p.adjust(Results$Pvalue,method = "fdr"))
Results <- cbind(Results,Region = rep(Region,dim(Results)[1]), Comparison = rep(paste(group1, group2, sep="/"),dim(Results)[1]))

#Append to Final Result Table
AllResults <- rbind(AllResults,Results)

#Write final Results Table
write.table(AllResults,"StatisticalResults_SubregionsvsSubregions_Region.csv",sep=";",row.names =F)


# ============== dHC vs vHC for selected subregions==============:
AllResults <-NULL
SubRegion <- "DG"
group1 <- "vHC"
group2 <- "dHC"
SubSet <- Proteins[Proteins$SubRegion == SubRegion & (Proteins$Region == group1 | Proteins$Region == group2),]

FilterMatrix <- spread(SubSet[,c("PG.ProteinAccessions","R.FileName","PG.Quantity")], R.FileName, PG.Quantity)
keep <- rowSums(is.na(FilterMatrix)) <= 4
ProteinAccesions <- FilterMatrix[keep,1]

Results <- NULL
for(i in 1:length(ProteinAccesions)){
  TestSet <- SubSet[SubSet$PG.ProteinAccessions == ProteinAccesions[i],]
  ANOVA <- aov(formula = PG.Quantity ~ Region + Error(AnimalNumber), TestSet)
  pvalue <- as.numeric(unlist(summary(ANOVA))["Error: Within.Pr(>F)1"])
  log2FC <- as.numeric(mean(TestSet[TestSet$Region == group1,]$PG.Quantity, na.rm = T) - mean(TestSet[TestSet$Region == group2,]$PG.Quantity, na.rm = T))
  Results <- rbind(Results,data.frame(paste(ProteinAccesions[i]),pvalue, log2FC))
}


names(Results) <- c("Accession","Pvalue","Log2FC")
assign <- match(Results$Accession,AccessionNameMapping$PG.ProteinAccessions)
Results <- cbind(Results, ShortName = AccessionNameMapping[assign,"PG.Genes"])
Results <- cbind(Results,fdr = p.adjust(Results$Pvalue,method = "fdr"))
Results <- cbind(Results,SubRegion = rep(SubRegion,dim(Results)[1]), Comparison = rep(paste(group1, group2, sep="/"),dim(Results)[1]))

AllResults <- rbind(AllResults,Results)

write.table(AllResults,"StatisticalResults_RegionsvsRegions_WithinSubRegions.csv",sep=";",row.names =F)


# ======================Region*SubRegion using ANOVA within animal ======================
SubSet <- Proteins

# Create Matrix to filter out low abundante Proteins, Create Target list for Testing
FilterMatrix <- spread(SubSet[,c("PG.ProteinAccessions","R.FileName","PG.Quantity")], R.FileName, PG.Quantity)
keep <- rowSums(is.na(FilterMatrix)) <= 4
ProteinAccesions <- FilterMatrix[keep,1]

# Performa ANOVA for all proteins
Results <- NULL
ResultsWideFormat <- NULL
for(i in 1:length(ProteinAccesions)){
  TestSet <- SubSet[SubSet$PG.ProteinAccessions == ProteinAccesions[i],]
  ANOVA <- aov(formula = PG.Quantity ~ Region * SubRegion + Error(AnimalNumber), TestSet)
  pvalue.Region <- as.numeric(unlist(summary(ANOVA))["Error: Within.Pr(>F)1"])
  pvalue.SubRegion <- as.numeric(unlist(summary(ANOVA))["Error: Within.Pr(>F)2"])
  pvalue.Interaction <- as.numeric(unlist(summary(ANOVA))["Error: Within.Pr(>F)3"])
  log2FC.Region <- as.numeric(mean(TestSet[TestSet$Region == "vHC",]$PG.Quantity, na.rm = T) - mean(TestSet[TestSet$Region == "dHC",]$PG.Quantity, na.rm = T))
  log2FC.CA1CA3 <- as.numeric(mean(TestSet[TestSet$SubRegion == "CA1",]$PG.Quantity, na.rm = T) - mean(TestSet[TestSet$SubRegion == "CA3",]$PG.Quantity, na.rm = T))
  log2FC.CA1DG <- as.numeric(mean(TestSet[TestSet$SubRegion == "CA1",]$PG.Quantity, na.rm = T) - mean(TestSet[TestSet$SubRegion == "DG",]$PG.Quantity, na.rm = T))
  log2FC.CA3DG <- as.numeric(mean(TestSet[TestSet$SubRegion == "CA3",]$PG.Quantity, na.rm = T) - mean(TestSet[TestSet$SubRegion == "DG",]$PG.Quantity, na.rm = T))
  Results <- rbind(Results,data.frame(paste(ProteinAccesions[i]),pvalue = pvalue.Region, Comparison = "vHC/dHC", Within = "all Subregions", Test = "Region", log2FC = log2FC.Region))
  Results <- rbind(Results,data.frame(paste(ProteinAccesions[i]),pvalue = pvalue.SubRegion, Comparison = "CA1/CA3", Within = "both Regions", Test = "Subregion",log2FC = log2FC.CA1CA3))
  Results <- rbind(Results,data.frame(paste(ProteinAccesions[i]),pvalue = pvalue.SubRegion, Comparison = "CA1/DG", Within = "both Regions", Test = "Subregion", log2FC = log2FC.CA1DG))
  Results <- rbind(Results,data.frame(paste(ProteinAccesions[i]),pvalue = pvalue.SubRegion, Comparison = "CA3/DG", Within = "both Regions", Test = "Subregion", log2FC = log2FC.CA3DG))
  Results <- rbind(Results,data.frame(paste(ProteinAccesions[i]),pvalue = pvalue.Interaction, Comparison = NA, Within = NA, Test = "Interaction", log2FC = NA))
  ResultsWideFormat <- rbind(ResultsWideFormat, data.frame(Accession = paste(ProteinAccesions[i]), pvalue.Region, pvalue.SubRegion, pvalue.Interaction, log2FC.Region, log2FC.CA1CA3, log2FC.CA1DG, log2FC.CA3DG))
}

# Create Results table and perform FDR testing
names(Results) <- c("Accession","Pvalue","Comparison","Within","Test","Log2FC")
assign <- match(Results$Accession,AccessionNameMapping$PG.ProteinAccessions)
Results <- cbind(Results, ShortName = AccessionNameMapping[assign,"PG.Genes"])
assign <- match(ResultsWideFormat$Accession,AccessionNameMapping$PG.ProteinAccessions)
ResultsWideFormat <- cbind(ShortName = AccessionNameMapping[assign,"PG.Genes"],ResultsWideFormat)
ResultsWideFormat <- cbind(ResultsWideFormat,fdr.Region = p.adjust(ResultsWideFormat$pvalue.Region,method = "fdr"),fdr.SubRegion = p.adjust(ResultsWideFormat$pvalue.SubRegion,method = "fdr"),fdr.Interaction = p.adjust(ResultsWideFormat$pvalue.Interaction,method = "fdr"))


library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}


ResultsWideFormat$meanSubRegionLogFC <- (abs(ResultsWideFormat$log2FC.CA1CA3) + abs(ResultsWideFormat$log2FC.CA1DG)+ abs(ResultsWideFormat$log2FC.CA3DG)) / 3

ggplot(ResultsWideFormat,aes(meanSubRegionLogFC,pvalue.SubRegion, color = p.adjust(ResultsWideFormat$pvalue.SubRegion,method = "fdr") < 0.05))+
  geom_point()+
  scale_colour_manual(values = c("black","green4","grey30","red","salmon")) + 
  scale_y_continuous(trans=reverselog_trans(10)) + 
  theme_bw()

ResultsWideFormat$FDR.SubRegion <- p.adjust(ResultsWideFormat$pvalue.SubRegion, method = "fdr")
ResultsWideFormat$FDR.Region <- p.adjust(ResultsWideFormat$pvalue.Region, method = "fdr")
ResultsWideFormat$FDR.Interaction <- p.adjust(ResultsWideFormat$pvalue.Interaction, method = "fdr")

FCcutoff <- 0
fdrCutoff <- 0.05
change <- NULL
for(i in 1:dim(ResultsWideFormat)[1]){
  if(ResultsWideFormat[i,"log2FC.Region"] > FCcutoff & ResultsWideFormat[i,"FDR.Region"] <= fdrCutoff){
    change<-append(change,"up-regulated 5% FDR")
  }
  else if(ResultsWideFormat[i,"log2FC.Region"] < -FCcutoff & ResultsWideFormat[i,"FDR.Region"] <= fdrCutoff){
    change<-append(change,"down-regulated 5% FDR")
  }
  else{
    change<-append(change,"not significant")
  }
}
ResultsWideFormat$change_Region <- change


ggplot(ResultsWideFormat,aes(log2FC.Region,pvalue.Region, color = change_Region))+
  geom_point()+
  scale_colour_manual(values = c("blue","grey30","red")) + 
  scale_y_continuous(trans=reverselog_trans(10)) + 
  theme_bw()

ggplot(ResultsWideFormat,aes(pvalue.Interaction)) + geom_histogram(bins = 25, color = "black", fill = "grey") + theme_bw()

assign <- match(ResultsWideFormat$Accession,Proteins$PG.ProteinAccessions)
ResultsWideFormat$GeneSymbol <- Proteins[assign,"PG.Genes"]

write.table(ResultsWideFormat,"SubregionxRegionANOVAResults.csv", sep = ";", row.names = F)



#===========================================================================================================================================
#=======================================================TWO GROUP COMPARISONS CVs vs Transcriptomics=====================================================
#===========================================================================================================================================


# ============== Swim vs Homecage within selected region and sub-region==============:
AllResults <-NULL
for(i in unique(Proteins$SubRegion)){
  for(j in unique(Proteins$Region)){
    
    print(paste(i, j , sep = " "))
    
    SubSet <- Proteins[Proteins$SubRegion == i & Proteins$Region == j,]
    group1 = "Swim"
    group2 = "Homecage"
    
    # Create Matrix to filter out low abundante Proteins, Create Target list for Testing
    FilterMatrix <- spread(SubSet[,c("PG.ProteinAccessions","R.FileName","PG.Quantity")], R.FileName, PG.Quantity)
    keep <- rowSums(is.na(FilterMatrix)) <= 4
    ProteinAccesions <- FilterMatrix[keep,1]
    
    # Performa ANOVA for all proteins
    Results <- NULL
    for(k in 1:length(ProteinAccesions)){
      TestSet <- SubSet[SubSet$PG.ProteinAccessions == ProteinAccesions[k],]
      CV = sd(TestSet$RawQuant) / mean(TestSet$RawQuant)
      Results <- rbind(Results,data.frame(Accession = paste(ProteinAccesions[k]), CV = CV, logQuant = mean(TestSet$PG.Quantity)))
    }
    
    # Create Results table and perform FDR testing
    assign <- match(Results$Accession,AccessionNameMapping$PG.ProteinAccessions)
    Results <- cbind(Results, ShortName = AccessionNameMapping[assign,"PG.Genes"])
    Results <- cbind(Results,Region = j)
    Results <- cbind(Results,SubRegion = i)
    
    #Append to Final Result Table
    AllResults <- rbind(AllResults, Results)
  }
}
ggplot(AllResults,aes(logQuant, CV)) + geom_point(size = 0.2) + stat_smooth() + geom_hline(yintercept = mean(CV), color = "red") + facet_grid(Region~SubRegion) + theme_bw()

GeneCounts <- read.table("P:/Lukas/Sequencing/20180813_TimeSeriesFST/GeneCountsTimeSeriesFST.txt")
GeneCounts$Group <- factor(GeneCounts$Group, levels = c("Homecage","Swim_45min","Swim_1h30min","Swim_2h","Swim_3h","Swim_4h"))
GeneCounts$Region <- factor(GeneCounts$Region, levels = c("vHC","dHC"))

s2c <- read.table(file.path("P:/Lukas/Sequencing/20180813_TimeSeriesFST/metadata", "TimeSeriesFST_s2c.txt"), header = TRUE, stringsAsFactors=FALSE, na.strings = "<NA>")


Regions <- c("vHC","dHC")
Groups <-c("Swim_45min","Swim_1h30min","Swim_2h","Swim_3h","Swim_4h")
CVs <- NULL

pdf("CVPlots.pdf", width = 6, height = 4)
for(i in unique(AllResults$Region)){
  for(j in unique(AllResults$SubRegion)){
    subset <- AllResults[(AllResults$Region == i & AllResults$SubRegion == j),]
    p1 <- ggplot(subset,aes(logQuant, CV)) + 
      geom_point(size = 0.5) + 
      stat_smooth() + 
      geom_hline(yintercept = mean(na.omit(subset$CV)), color = "red") + 
      theme_bw() + 
      ggtitle(paste("Proteomics",i,j, sep = " ")) + 
      ylab("Biological coefficient of variation") + xlab("log Intensity")
    print(p1)
    CVs <- rbind(CVs, data.frame(meanCV = mean(na.omit(subset$CV)), method = "Proteomics"))
  }
}

for(i in Regions){
    GeneCountsSubset <-  GeneCounts[GeneCounts$Region == i,]
    
    GeneCountsSubset <- data.frame(ensembl_gene_id = GeneCountsSubset$ensembl_gene_id,sample = GeneCountsSubset$sample ,counts = GeneCountsSubset$est_counts)
    GeneCountsSubset <- spread(GeneCountsSubset, sample, counts)
    row.names(GeneCountsSubset) <- GeneCountsSubset[,1]
    Conditions <- s2c[match(colnames(GeneCountsSubset)[2:ncol(GeneCountsSubset)],s2c$SampleName), "Condition",]
    
    
    y <- DGEList(counts=GeneCountsSubset[,c(2:ncol(GeneCountsSubset))], group = Conditions)
    y <- calcNormFactors(y)
    y <- estimateDisp(y)
    plotBCV(y, main=paste("Transcriptomics", i,j, sep = " "))
    CVs <- rbind(CVs, data.frame(meanCV = sqrt(y$common.disp), method = "Transcriptomics"))
  }

p1 <- ggplot(CVs, aes(method,meanCV)) + geom_boxplot() + geom_point() + theme_bw()
print(p1)

power <- NULL

ProtCV <- mean(CVs$meanCV[CVs$method == "Proteomics"])
TransCV <- mean(CVs$meanCV[CVs$method == "Transcriptomics"])

for(delta in seq(from = 0.05, to = 1, by = 0.05)){
  power.t.test(n = 6,sd = 0.1375,sig.level = 0.05, delta = delta)
  power <- rbind(power, data.frame(power = power.t.test(n = 8,sd = ProtCV,sig.level = 0.05, delta = delta)$power, delta = delta, type = "Proteomics"))
  power <- rbind(power, data.frame(power = power.t.test(n = 7.5,sd = TransCV,sig.level = 0.05, delta = delta)$power, delta = delta, type = "Transcriptomics"))
}
p1 <- ggplot(power,aes(delta,power,color = type)) + geom_line() + theme_bw()
print(p1)

dev.off()

