.libPaths(append(.libPaths(),"C:/Users/lukasv/Documents/RUserLibs"))

#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR","C:/Users/lukasv/Documents/RUserLibs")
#install.packages("aggregation","C:/Users/lukasv/Documents/RUserLibs")


library(ImpulseDE2)
library(ggplot2)
library(edgeR)
library(tidyr)
library(aggregation)


setwd("P:/Lukas/Sequencing/20180813_TimeSeriesFST")


# Set WorkinDir and read Sample to Conditions mapping table
sample_id <- dir(file.path("kallisto"))
s2c <- read.table(file.path("metadata", "TimeSeriesFST_s2c.txt"), header = TRUE, stringsAsFactors=FALSE, na.strings = "<NA>")
s2c <- dplyr::mutate(s2c, path = file.path("kallisto", sample_id, "output"))

# Create IDsmapping list
IDsmapping <- read.table(file = paste(s2c$path[1],"/abundance.tsv",sep = ''), sep = '\t', header = TRUE)
IDsmapping <- data.frame(IDsmapping$target_id, do.call('rbind', strsplit(as.character(IDsmapping$target_id),'|',fixed=TRUE)))
names(IDsmapping) <-c("ID","ensembl_transcript_id","ensembl_gene_id","VEGA_gene_id","VEGA_transcript_id","ensembl_transcript_name","ensembl_gene_name","length","type")
IDsmapping$ensembl_gene_id <- substring(IDsmapping$ensembl_gene_id,1,18)

TranscriptCounts <- NULL
GeneCounts <- NULL

for(i in 1:dim(s2c)[1]){
  print(paste(i,s2c[i,1],sep = ": "))
  x <- read.table(file = paste(s2c$path[i],"/abundance.tsv",sep = ''), sep = '\t', header = TRUE)
  sample <- rep(s2c[i,1],dim(x)[1])
  x <- cbind(x,sample)
  
  #Append to TranscriptCounts list
  assign <- match(x$target_id,IDsmapping$ID)
  x <- data.frame(ensembl_transcript_id = substring(IDsmapping$ensembl_transcript_id[assign],1,18), 
                  ensembl_gene_id = substring(IDsmapping$ensembl_gene_id[assign],1,18), 
                  ensembl_transcript_name = IDsmapping$ensembl_transcript_name[assign], 
                  ensembl_gene_name =IDsmapping$ensembl_gene_name[assign], 
                  type = IDsmapping$type[assign],
                  x[,2:dim(x)[2]])
  
  assign <- match(x$sample,s2c$SampleName)
  x <- cbind(x, SampleProcessing = s2c[assign,"SampleProcessing"],
             CircadianTime = s2c[assign,"CircadianTime"], 
             Region = s2c[assign,"Region"], 
             Condition = s2c[assign,"Condition"], 
             Group = s2c[assign,"Group"], 
             CageNumber = s2c[assign,"CageNumber"], 
             Block = s2c[assign,"Block"])
  
  #Aggregate on the gene level and append to GeneCounts list
  keep <- x$type == "protein_coding"
  y <- data.frame(ensembl_gene_id = x$ensembl_gene_id,est_counts = x$est_counts,tpm = x$tpm)
  y <- y[keep,]
  y <- aggregate(y[,c(2,3)], by = list(unique.values = y$ensembl_gene_id) ,FUN = sum)
  
  assign <- match(y$unique.values,x$ensembl_gene_id)
  y <- data.frame(ensembl_gene_id = x$ensembl_gene_id[assign], 
                  ensembl_gene_name =x$ensembl_gene_name[assign],
                  y[,c(2,3)],
                  sample = x$sample[assign],
                  SampleProcessing = x$SampleProcessing[assign],
                  CircadianTime = x$CircadianTime[assign],
                  Region = x$Region[assign],
                  Condition = x$Condition[assign],
                  Group = x$Group[assign],
                  CageNumber = x$CageNumber[assign],
                  Block = x$Block[assign])
  
  TranscriptCounts <- rbind(TranscriptCounts,x)
  GeneCounts <- rbind(GeneCounts,y)
}


#Explore Transcripts in the time series
TranscriptCounts <- read.table("TranscriptCountsTimeSeriesFST.txt")
TranscriptCounts$Group <- factor(TranscriptCounts$Group, levels = c("Homecage","Swim_45min","Swim_1h30min","Swim_2h","Swim_3h","Swim_4h"))
TranscriptCounts$Region <- factor(TranscriptCounts$Region, levels = c("vHC","dHC"))
GeneCounts <- read.table("GeneCountsTimeSeriesFST.txt")
GeneCounts$Group <- factor(GeneCounts$Group, levels = c("Homecage","Swim_45min","Swim_1h30min","Swim_2h","Swim_3h","Swim_4h"))
GeneCounts$Region <- factor(GeneCounts$Region, levels = c("vHC","dHC"))



TargetList <- rbind("Cyr61")
Targets <- NULL
for(i in TargetList){
  Targets <- rbind(Targets,GeneCounts[GeneCounts$ensembl_gene_name == i,])
}
#Targets <- rbind(Targets[Targets$ensembl_transcript_id =="ENSMUST00000115567",],Targets[Targets$ensembl_transcript_id =="ENSMUST00000115571",])

ggplot(Targets,aes(ensembl_gene_name,tpm, color=Group)) + 
  geom_boxplot(outlier.shape = NA) + 
  facet_grid(Region~.)+  
  geom_point(position = position_jitterdodge())+ 
  theme_bw() + 
  ggtitle(paste("Expression of ",Targets$ensembl_gene_name,sep ="")) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values = c("grey28","#CE2576","#D05172","#D37366","#D68C4F","#DAA521"))

#EdgeR analysis for individual time points

Regions <- c("vHC","dHC")
Groups <-c("Swim_45min","Swim_1h30min","Swim_2h","Swim_3h","Swim_4h")

AllResults <- NULL
for(i in Regions){
  for(j in Groups){
    GeneCountsSubset <-  GeneCounts[GeneCounts$Region == i & (GeneCounts$Group == j | GeneCounts$Group =="Homecage"),]
    
    GeneCountsSubset <- data.frame(ensembl_gene_id = GeneCountsSubset$ensembl_gene_id,sample = GeneCountsSubset$sample ,counts = GeneCountsSubset$est_counts)
    GeneCountsSubset <- spread(GeneCountsSubset, sample, counts)
    row.names(GeneCountsSubset) <- GeneCountsSubset[,1]
    Conditions <- s2c[match(colnames(GeneCountsSubset)[2:ncol(GeneCountsSubset)],s2c$SampleName), "Condition",]
    
    
    y <- DGEList(counts=GeneCountsSubset[,c(2:ncol(GeneCountsSubset))], group = Conditions)
    y <- calcNormFactors(y)
    y <- estimateDisp(y)
    
    keep <- rowSums(y$counts>10) >= 5
    y <- y[keep, , keep.lib.sizes=FALSE]
    
    Results <- exactTest(y)$table
    Results <- cbind(Results,p.adjust(Results[,3],method = "fdr"))
    assign <- match(row.names(Results),IDsmapping$ensembl_gene_id)
    Results <- cbind(Results, IDsmapping[assign,"ensembl_gene_id"], IDsmapping[assign,"ensembl_gene_name"])
    names(Results) <- c("logFC","logCPM","P.Value","adj.P.Val","EnsembleID","ShortName")
    AllResults <- rbind(AllResults, cbind(Results, Region = i, TimePoint = j, Test= paste(i,j,sep="_")))
  }
}
write.table(AllResults,"IndividualTimepoints.csv", sep = ";",row.names = F)
Plot.Volcanos.MultiFDR(AllResults,path = paste(getwd(),"IndividualTimePoints_Volcanoplots.pdf",sep ="/"), width = 12, height = 6)


#EdgeR analysis over timeseries in vHC and dHC independently

Regions <- c("vHC","dHC")

AllResults <- NULL
for(i in Regions){

GeneCountsSubset <- GeneCounts[GeneCounts$Region == i,]
GeneCountsSubset <- data.frame(ensembl_gene_id = GeneCountsSubset$ensembl_gene_id,sample = GeneCountsSubset$sample ,counts = GeneCountsSubset$est_counts)
GeneCountsSubset <- spread(GeneCountsSubset, sample, counts)
row.names(GeneCountsSubset) <- GeneCountsSubset[,1]
Groups = s2c[match(colnames(GeneCountsSubset)[2:ncol(GeneCountsSubset)],s2c$SampleName), "Group",]
design <- model.matrix(~Groups)

y <- DGEList(counts=GeneCountsSubset[,c(2:ncol(GeneCountsSubset))], group = Groups)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)

keep <- rowSums(y$counts>10) >= 5
y <- y[keep, , keep.lib.sizes=FALSE]

fit <- glmQLFit(y,design)
Results <- glmQLFTest(fit, colnames(design)[-1])
Results <- Results$table
Results$FDR <- p.adjust(Results[,"PValue"],method = "fdr")
assign <- match(row.names(Results),IDsmapping$ensembl_gene_id)
Results <- cbind(Results, EnsembleID = IDsmapping[assign,"ensembl_gene_id"], ShortName = IDsmapping[assign,"ensembl_gene_name"])
Results$Region <- i
AllResults <- rbind(Results,AllResults)
}

write.table(Results,"dHC_Groups.csv",sep =";",row.names = F)


ggplot(AllResults,aes(PValue))+
  geom_histogram(color = "black", fill = "grey90",breaks = seq(from = 0, to = 1, by = 0.01)) +
  theme_bw()+
  facet_grid(Region~.)


# Create log2FC expression matrixes 
vHC <- read.table("vHC_IndividualTimepoints.csv",sep=";")
dHC <- read.table("dHC_IndividualTimepoints.csv",sep=";")

vHCPs <- spread(vHC[,c("EnsembleID","Pvalue","Test")], Test, Pvalue)
vHCPs <- vHCPs[,c("EnsembleID","Swim_45min","Swim_1h30min","Swim_2h","Swim_3h","Swim_4h")]
dHCPs <- spread(dHC[,c("EnsembleID","Pvalue","Test")], Test, Pvalue)
dHCPs <- dHCPs[,c("EnsembleID","Swim_45min","Swim_1h30min","Swim_2h","Swim_3h","Swim_4h")]


vHCLogFCs <- spread(vHC[,c("EnsembleID","logFC","Test")], Test, logFC)
vHCLogFCs <- vHCLogFCs[,c("EnsembleID","Swim_45min","Swim_1h30min","Swim_2h","Swim_3h","Swim_4h")]
dHCLogFCs <- spread(dHC[,c("EnsembleID","logFC","Test")], Test, logFC)
dHCLogFCs <- dHCLogFCs[,c("EnsembleID","Swim_45min","Swim_1h30min","Swim_2h","Swim_3h","Swim_4h")]


#Pvalue aggregation 
vHCaggregatedPs <- NULL
for(i in 1:dim(vHCPs)[1]){
  vHCaggregatedPs <- rbind(vHCaggregatedPs,data.frame(EnsembleGeneID = vHCPs[i,1],AggregatedPValue = fisher(vHCPs[i,2:5])))
}
dHCaggregatedPs <- NULL
for(i in 1:dim(dHCPs)[1]){
  dHCaggregatedPs <- rbind(dHCaggregatedPs,data.frame(EnsembleGeneID = dHCPs[i,1],AggregatedPValue = fisher(dHCPs[i,2:5])))
}

vHCaggregatedPs <- cbind(vHCaggregatedPs, FDR = p.adjust(vHCaggregatedPs$AggregatedPValue,method = "fdr"))
vHCaggregatedPs <- vHCaggregatedPs[(vHCaggregatedPs$FDR <= 0.05),]
assign <- match(vHCaggregatedPs$EnsembleGeneID,vHC$EnsembleID)
vHCaggregatedPs <- cbind(ShortName = vHC[assign,"ShortName"], vHCaggregatedPs)
dHCaggregatedPs <- cbind(dHCaggregatedPs, FDR = p.adjust(dHCaggregatedPs$AggregatedPValue,method = "fdr"))
dHCaggregatedPSig <- dHCaggregatedPs[(dHCaggregatedPs$FDR <= 0.05),]
assign <- match(dHCaggregatedPs$EnsembleGeneID,dHC$EnsembleID)
dHCaggregatedPs <- cbind(ShortName = dHC[assign,"ShortName"], dHCaggregatedPs)

write.table(vHCaggregatedPs,"vHCaggregatedFDR.csv",sep =";",row.names = F)
write.table(dHCaggregatedPs,"dHCaggregatedFDR.csv",sep =";",row.names = F)



assign <- match(vHCaggregatedPs$EnsembleGeneID,vHCLogFCs$EnsembleID)
vHCSigLogFCs <- vHCLogFCs[assign,]
write.table(na.omit(vHCSigLogFCs),"vHCLogFCmatrix_Sig.txt",sep ="\t",row.names = F)
write.table(na.omit(dHCLogFCs),"dHCLogFCmatrix.txt",sep ="\t",row.names = F)


#create FDR matrixes

vHCFDRs <- spread(vHC[,c("EnsembleID","FDR","Test")], Test, FDR)
vHCFDRs <- vHCFDRs[,c("EnsembleID","Swim_45min","Swim_1h30min","Swim_2h","Swim_3h","Swim_4h")]
dHCFDRs <- spread(dHC[,c("EnsembleID","FDR","Test")], Test, FDR)
dHCFDRs <- dHCFDRs[,c("EnsembleID","Swim_45min","Swim_1h30min","Swim_2h","Swim_3h","Swim_4h")]

assign <- match(vHCFDRs$EnsembleID,vHC$EnsembleID)
vHCFDRs <- cbind(ShortName = vHC[assign,"ShortName"],vHCFDRs)
assign <- match(dHCFDRs$EnsembleID,dHC$EnsembleID)
dHCFDRs <- cbind(ShortName = dHC[assign,"ShortName"],dHCFDRs)

write.table(vHCFDRs,"vHCFDRmatrix.csv",sep =";",row.names = F)
write.table(dHCFDRs,"dHCFDRmatrix.csv",sep =";",row.names = F)

