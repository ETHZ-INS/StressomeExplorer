.libPaths(append(.libPaths(),"C:/Users/lukasv/Documents/RUserLibs"))

setwd("P:/Lukas/Analysis/20191104_AllPhosphodata/")

library(tidyr)
library(lme4)
library(lmerTest)
library(DEP)
library(ggplot2)
library(curl)
library(aggregation)
library(stringr)
library(cowplot)
library(SummarizedExperiment)


data <- read.table("20191031_174447_20191028_AllPhosphoData_Report.csv",sep = ",", header = T)
s2c <- read.table("metadata/s2c_allphos.csv", sep = ";", header = T)
rownames(s2c) <- s2c$FileName

length(unique(data$EG.PrecursorId))

matr <- data[,c("R.FileName","EG.ModifiedPeptide","FG.Quantity")]
matr <- aggregate(FG.Quantity ~ R.FileName + EG.ModifiedPeptide, data = matr, sum)


matr <- spread(matr, key = R.FileName, value = FG.Quantity)

rownames(matr) <- matr[,1]
matr <- as.matrix(matr[,-1])
matr <- log2(matr)

#homecage vs immediate within regions
AllResults <- NULL
minpergroup <- 3
mintotal <- 10

for(i in unique(s2c$Region)){
  submat <- matr[,s2c[colnames(matr),"Region"] == i & s2c[colnames(matr),"Experiment"] == "6minFST"]
  Condition <- s2c[colnames(submat),"Condition"]
  groups <- model.matrix(~Condition)[,2]
  keep <- apply(submat,1,function(x){
    min(sum(!is.na(x[as.logical(groups)])),sum(!is.na(x[!as.logical(groups)])))
  })
  submat <- submat[keep,]
  
  Results <- NULL
  for(j in 1:nrow(submat)){
    tryCatch({
    testdat <- submat[j,]
    model <- lm(groups~testdat)
    FC <- mean(testdat[as.logical(groups)], na.rm = T) - mean(testdat[!as.logical(groups)], na.rm = T)
    Results <- rbind(Results, data.frame(EG = rownames(submat)[j], p.value = summary(model)$coefficients[,"Pr(>|t|)"][2], logFC = FC, Region = i))
    }, error = function(e){print(e)})
  }
  AllResults <- rbind(AllResults,data.frame(Results, FDR = p.adjust(Results$p.value, method = "fdr")))
}

write.table(AllResults,"ImmediatevsControl.csv", sep = ";", row.names = F)


#homecage vs timepoints within regions
minpergroup <- 3

for(i in unique(s2c$Region)){
  for(k in c("15min","30min","45min")){
  submat <- matr[,s2c[colnames(matr),"Region"] == i & s2c[colnames(matr),"Experiment"] == "TimeSeriesFST" & s2c[colnames(matr),"TimePoint"] %in% c(k,"0min")]
  Condition <- s2c[colnames(submat),"Condition"]
  groups <- model.matrix(~Condition)[,2]
  keep <- apply(submat,1,function(x){
    min(sum(!is.na(x[as.logical(groups)])),sum(!is.na(x[!as.logical(groups)]))) >= minpergroup
  })
  submat <- submat[keep,]
  
  Results <- NULL
  for(j in 1:nrow(submat)){
    tryCatch({
      testdat <- submat[j,]
      model <- lm(groups~testdat)
      FC <- mean(testdat[as.logical(groups)], na.rm = T) - mean(testdat[!as.logical(groups)], na.rm = T)
      Results <- rbind(Results, data.frame(EG = rownames(submat)[j], p.value = summary(model)$coefficients[,"Pr(>|t|)"][2], logFC = FC, Region = i, TimePoint = k))
    }, error = function(e){print(e)})
  }
  AllResults <- rbind(AllResults,data.frame(Results, FDR = p.adjust(Results$p.value, method = "fdr")))
}
}



write.table(AllResults,"AllComps_LM.csv", sep = ";", row.names = F)
AllResults$Condition2 <- apply(AllResults,1,FUN = function(x){paste(x[4],x[6], sep = "_")})
Resmat <- spread(AllResults[,c("EG","logFC","Condition2")], key = Condition2, value = logFC)


#LMER with blocking immediate vs control
AllResults <- NULL
minpergroup <- 3
mintotal <- 10

for(i in unique(s2c$Region)){
  submat <- matr[,s2c[colnames(matr),"Region"] == i & s2c[colnames(matr),"Experiment"] == "6minFST"]
  Condition <- s2c[colnames(submat),"Condition"]
  Replicate <- as.factor(s2c[colnames(submat),"Replicate"])
  groups <- model.matrix(~Condition)[,2]
  keep <- apply(submat,1,function(x){
    min(sum(!is.na(x[as.logical(groups)])),sum(!is.na(x[!as.logical(groups)]))) >= minpergroup & sum(!is.na(x)) >= mintotal
  })
  submat <- submat[keep,]
  
  Results <- NULL
  for(j in 1:nrow(submat)){
    if(j %% 1000 == 0){
      print(paste(i,j, sep =" "))
    }
    tryCatch({
      testdat <- data.frame(logIntensity = submat[j,], groups = groups, Replicate = Replicate)
      model <- lmer(data = testdat, formula = logIntensity~groups + (1|Replicate))
      
      FC <- mean(testdat[as.logical(groups),"logIntensity"], na.rm = T) - mean(testdat[!as.logical(groups),"logIntensity"], na.rm = T)
      Results <- rbind(Results, data.frame(EG = rownames(submat)[j], p.value = summary(model)$coefficients[,"Pr(>|t|)"][2], logFC = FC, Region = i))
    }, error = function(e){print(e)})
  }
  AllResults <- rbind(AllResults,data.frame(Results, FDR = p.adjust(Results$p.value, method = "fdr")))
}

write.table(AllResults,"ImmediatevsControl_Lmer.csv", sep = ";", row.names = F)

#LMER with blocking timepoints vs control


AllResults <- NULL
minpergroup <- 3
mintotal <- 7

for(i in unique(s2c$Region)){
  for(k in c("15min","30min","45min")){
    submat <- matr[,s2c[colnames(matr),"Region"] == i & s2c[colnames(matr),"Experiment"] == "TimeSeriesFST" & s2c[colnames(matr),"TimePoint"] %in% c(k,"0min")]
    Condition <- s2c[colnames(submat),"Condition"]
    Replicate <- as.factor(s2c[colnames(submat),"Replicate"])
    groups <- model.matrix(~Condition)[,2]
    keep <- apply(submat,1,function(x){
      min(sum(!is.na(x[as.logical(groups)])),sum(!is.na(x[!as.logical(groups)]))) >= minpergroup & sum(!is.na(x)) >= mintotal
    })
    submat <- submat[keep,]
    
    Results <- NULL
  for(j in 1:nrow(submat)){
    if(j %% 1000 == 0){
      print(paste(i,k,j, sep =" "))
    }
    tryCatch({
      testdat <- data.frame(logIntensity = submat[j,], groups = groups, Replicate = Replicate)
      model <- lmer(data = testdat, formula = logIntensity~groups + (1|Replicate))
      
      FC <- mean(testdat[as.logical(groups),"logIntensity"], na.rm = T) - mean(testdat[!as.logical(groups),"logIntensity"], na.rm = T)
      Results <- rbind(Results, data.frame(EG = rownames(submat)[j], p.value = summary(model)$coefficients[,"Pr(>|t|)"][2], logFC = FC, Region = i, TimePoint = k))
    }, error = function(e){print(e)})
  }
  }
  AllResults <- rbind(AllResults,data.frame(Results, FDR = p.adjust(Results$p.value, method = "fdr")))
}

write.table(AllResults,"TimePointsvsControl_Lmer.csv", sep = ";", row.names = F)


#Run DEP for combined data

SE <- SummarizedExperiment( list(),
                            colData=s2c,
                            rowDat=)



data <- read.table("20191031_174447_20191028_AllPhosphoData_Report.csv",sep = ",", header = T)
s2c <- read.table("metadata/s2c_allphos.csv", sep = ";", header = T)
rownames(s2c) <- s2c$FileName

length(unique(data$EG.PrecursorId))

matr <- data[,c("R.FileName","EG.ModifiedPeptide","FG.Quantity")]
matr <- aggregate(FG.Quantity ~ R.FileName + EG.ModifiedPeptide, data = matr, sum)


matr <- spread(matr, key = R.FileName, value = FG.Quantity)

rownames(matr) <- matr[,1]
matr <- matr[,-1]


s2csub <- s2c
mat <- matr[,names(matr) %in% s2csub$FileName]
mat$Protein.names <- row.names(mat)
mat$Protein.IDs <- paste("ID_",c(1:nrow(mat)),sep = "")
data_unique <- make_unique(mat, "Protein.names", "Protein.IDs", delim = ";")

experimental_design <- data.frame(label = as.character(s2csub[,"FileName"]), 
                                  condition = as.character(paste(s2csub[,"ShortExp"],s2csub[,"Region"],s2csub[,"ShortCond"], sep = "_")), 
                                  replicate = as.integer(s2csub[,"Replicate"]))
experimental_design[,c("label","condition")] <- lapply(experimental_design[,c("label","condition")], as.character)


#Create DEP data
LFQ_columns <- na.omit(match(s2csub$FileName,colnames(data_unique)))
data_se <- make_se(data_unique, LFQ_columns, experimental_design)
data_filt <- filter_missval(data_se, thr = 2)
data_norm <- normalize_vsn(data_filt)
assays(data_filt)$vsn <- assay(data_norm)
data_imp <- impute(data_filt, fun = "MinProb", q = 0.01)
assays(data_filt)$imputed <- assay(data_imp)

data_diff_manual <- test_diff(data_imp, type = "manual", test = c("Exp1_vHC_SW06_vs_Exp1_vHC_CON",
                                                                  "Exp2_vHC_SW15_vs_Exp2_vHC_CON",
                                                                  "Exp2_vHC_SW30_vs_Exp2_vHC_CON",
                                                                  "Exp2_vHC_SW45_vs_Exp2_vHC_CON",
                                                                  "Exp1_dHC_SW06_vs_Exp1_dHC_CON",
                                                                  "Exp2_dHC_SW15_vs_Exp2_dHC_CON", 
                                                                  "Exp2_dHC_SW30_vs_Exp2_dHC_CON",
                                                                  "Exp2_dHC_SW45_vs_Exp2_dHC_CON"))


dep_man <- add_rejections(data_diff_manual, alpha = 0.05, lfc = log2(0.5))


#Export and build Results data.frame
data_results <- get_results(dep_man)
assign <- match(data_results$name,data$EG.ModifiedPeptide)
data_results$Protein <- data[assign,"PG.ProteinAccessions"]
data_results$ModifiedSequence <- data[assign,"EG.PTMLocalizationProbabilities"]
data_results$StrippedSequence <- data[assign,"EG.StrippedSequence"]
data_results$Positions <- data[assign,"EG.PTMPositions..Phospho_STY_NL_98..STY.."]
data_results$Probabilities <- data[assign,"EG.PTMProbabilities..Phospho_STY_NL_98..STY.."]

#Produce PDF with DEP plots

pdf("DEPPlots_AllComparisons.pdf", width = 20, height = 15)
plot_frequency(data_se)
plot_numbers(data_filt)
plot_coverage(data_filt)
plot_normalization(data_filt, data_norm)
plot_detect(data_filt)
plot_imputation(data_norm, data_imp)
plot_cor(dep_man, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
plot_heatmap(dep_man, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))  
plot_heatmap(dep_man, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)
plot_pca(dep_man, x = 1, y = 2, n = 1000, point_size = 4)

dev.off()

#Build up SE object
RD <- rowData(data_filt)
assign <- match(RD$name,data$EG.ModifiedPeptide)
RD$Protein <- data[assign,"PG.ProteinAccessions"]
RD$StrippedSequence <- data[assign,"EG.StrippedSequence"]
RD$Positions <- data[assign,"EG.PTMPositions..Phospho_STY_NL_98..STY.."]
RD$Probabilities <- data[assign,"EG.PTMProbabilities..Phospho_STY_NL_98..STY.."]

Mapping <- read.table("MappingTable.csv", header = T, sep = ";")
assign <- match(RD$Protein,Mapping$UniprotIDentifier)
RD$GeneSymbol <- Mapping[assign,"GeneName"]

assign <- match(RD$StrippedSequence, peppos$StripedSequence)
RD$EarliestPos <- peppos[assign,"earliestPos"]

PTMSites <- NULL
for(i in 1:nrow(RD)){
  if(RD[i,"Probabilities"] == ""){
    if(str_count(RD[i,"name"], pattern = "9") == 0){
      PTMSites <- append(PTMSites, "None")
    }else{
      PTMSites <- append(PTMSites, "Ambiguous")
    }
  }else{
    PTMSites <- append(PTMSites, getpositions(RD[i,"Positions"],RD[i,"EarliestPos"]))
  }
}
RD$PhosphoSites <- PTMSites
RD$NPhospho <- str_count(RD[,"name"], pattern = "9")

peppos <- NULL
for(i in unique(RD$Protein)){
  print(i)
  peppos <- rbind(peppos,pepAlign(unique(RD[RD$Protein == i, "StrippedSequence"]),str_split(i,pattern = ";")[[1]][1]))
}
RD$Protein.names <- NULL
RD$Peptide.ID <- RD$Protein.IDs
RD$Protein.IDs <- NULL
rowData(data_filt) <- RD

CD <- colData(data_filt)
assign <- match(CD$label,s2c$FileName)
CD$condition <- NULL
CD$Region <- s2c[assign,"Region"]
CD$Condition2 <- s2c[assign,"Condition2"]
CD$Animal <- s2c[assign,"Animal"]
CD$TimePoint <- s2c[assign,"TimePoint"]
CD$Block <- s2c[assign,"Block"]
CD$Condition <- s2c[assign,"Condition"]

CD$Experiment <- s2c[assign,"Experiment"]
colData(data_filt) <- CD

saveRDS(data_filt, "PhosphoData.SE.rds")


#Get positions and add phosphosites for results

pepAlign <- function(pep_seq, protein_id){
  prot <- paste(readLines(curl(paste0("https://www.uniprot.org/uniprot/",protein_id,".fasta")))[-1],collapse="")
  positions <- str_locate_all(prot, as.character(pep_seq))
  res <- data.frame( StripedSequence=pep_seq, 
                     pos_status=sapply(positions, FUN=function(x){ 
                       if(nrow(x)==0) return("not found")
                       ifelse(nrow(x)>1,"multimap","unique")
                     }),
                     earliestPos=sapply(positions, FUN=function(x){ if(nrow(x)==0) return(NA); min(x[,1]) }),
                     allPos=sapply(positions, FUN=function(x){ paste(paste(x[,1],x[,2],sep="-"), collapse="; ") })
  )
  res[order(res$earliestPos),]
}


peppos <- NULL
for(i in unique(data_results$Protein)){
  print(i)
  peppos <- rbind(peppos,pepAlign(unique(data_results[data_results$Protein == i, "StrippedSequence"]),str_split(i,pattern = ";")[[1]][1]))
}

assign <- match(data_results$StrippedSequence, peppos$StripedSequence)
data_results$EarliestPos <- peppos[assign,"earliestPos"]

getpositions <- function(x, ep){
  splitted <- str_split(x, pattern = ";")
  pos <- as.numeric(splitted[[1]]) + (ep - 1)
  paste(pos,sep = ";", collapse = ";")
}

PTMSites <- NULL
for(i in 1:nrow(data_results)){
  if(data_results[i,"Probabilities"] == ""){
    PTMSites <- append(PTMSites, "None")
  }else{
    PTMSites <- append(PTMSites, getpositions(data_results[i,"Positions"],data_results[i,"EarliestPos"]))
  }
}
data_results$PhosphoSites <- PTMSites

Mapping <- read.table("MappingTable.csv", header = T, sep = ";")
assign <- match(data_results$Protein,Mapping$UniprotIDentifier)
data_results$GeneSymbol <- Mapping[assign,"GeneName"]

data_results$Exp1_dHC_SW06_vs_Exp1_dHC_CON_FDR <- p.adjust(data_results$Exp1_dHC_SW06_vs_Exp1_dHC_CON_p.val, method = "fdr")
data_results$Exp1_vHC_SW06_vs_Exp1_vHC_CON_FDR <- p.adjust(data_results$Exp1_vHC_SW06_vs_Exp1_vHC_CON_p.val, method = "fdr")
data_results$Exp2_dHC_SW15_vs_Exp2_dHC_CON_FDR <- p.adjust(data_results$Exp2_dHC_SW15_vs_Exp2_dHC_CON_p.val, method = "fdr")
data_results$Exp2_dHC_SW30_vs_Exp2_dHC_CON_FDR <- p.adjust(data_results$Exp2_dHC_SW30_vs_Exp2_dHC_CON_p.val, method = "fdr")
data_results$Exp2_dHC_SW45_vs_Exp2_dHC_CON_FDR <- p.adjust(data_results$Exp2_dHC_SW45_vs_Exp2_dHC_CON_p.val, method = "fdr")
data_results$Exp2_vHC_SW15_vs_Exp2_vHC_CON_FDR <- p.adjust(data_results$Exp2_vHC_SW15_vs_Exp2_vHC_CON_p.val, method = "fdr")
data_results$Exp2_vHC_SW30_vs_Exp2_vHC_CON_FDR <- p.adjust(data_results$Exp2_vHC_SW30_vs_Exp2_vHC_CON_p.val, method = "fdr")
data_results$Exp2_vHC_SW45_vs_Exp2_vHC_CON_FDR <- p.adjust(data_results$Exp2_vHC_SW45_vs_Exp2_vHC_CON_p.val, method = "fdr")
data_results$AnySignificant <- apply(data_results[,c("Exp1_dHC_SW06_vs_Exp1_dHC_CON_FDR","Exp1_vHC_SW06_vs_Exp1_vHC_CON_FDR","Exp2_dHC_SW15_vs_Exp2_dHC_CON_FDR","Exp2_dHC_SW30_vs_Exp2_dHC_CON_FDR","Exp2_dHC_SW45_vs_Exp2_dHC_CON_FDR","Exp2_vHC_SW15_vs_Exp2_vHC_CON_FDR","Exp2_vHC_SW30_vs_Exp2_vHC_CON_FDR","Exp2_vHC_SW45_vs_Exp2_vHC_CON_FDR")],
      1,FUN = function(x){min(x) < 0.05})

#write.table(data_results,"DEPResults_AllComparisons.csv", sep = ";", row.names = F)

data_sig <- data_results[data_results$AnySignificant,]

data_sig6 <- data_sig
data_sig6$significance <- apply(data_sig6[,c("Exp1_dHC_SW06_vs_Exp1_dHC_CON_FDR","Exp1_vHC_SW06_vs_Exp1_vHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sig6 <- data_sig6[order(-data_sig6$significance),]
p1 <- ggplot(data_sig6,aes(Exp1_dHC_SW06_vs_Exp1_dHC_CON_ratio,Exp1_vHC_SW06_vs_Exp1_vHC_CON_ratio, color = -log10(data_sig6$significance))) + 
  geom_point() + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold"),
                         values = c(0, -log10(0.05)/max(-log10(data_sig6$significance)), 1)) +
  ggtitle("6 min Swim vs Control")+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(name = "dHC logFC Swim vs. Control") +
  scale_y_continuous(name = "vHC logFC Swim vs. Control") +
  labs(colour="-log10(FDR)") +
  theme_bw()

data_sig15 <- data_sig
data_sig15$significance <- apply(data_sig15[,c("Exp2_dHC_SW15_vs_Exp2_dHC_CON_FDR","Exp2_vHC_SW15_vs_Exp2_vHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sig15 <- data_sig15[order(-data_sig15$significance),]
p2 <- ggplot(data_sig15,aes(Exp2_dHC_SW15_vs_Exp2_dHC_CON_ratio,Exp2_vHC_SW15_vs_Exp2_vHC_CON_ratio, color = -log10(data_sig15$significance))) + 
  geom_point() + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/max(-log10(data_sig15$significance)), 1)) +
  ggtitle("15 min Swim vs Control")+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(name = "dHC logFC Swim vs. Control") +
  scale_y_continuous(name = "vHC logFC Swim vs. Control") +
  labs(colour="-log10(FDR)") +
  theme_bw()

data_sig30 <- data_sig
data_sig30$significance <- apply(data_sig30[,c("Exp2_dHC_SW30_vs_Exp2_dHC_CON_FDR","Exp2_vHC_SW30_vs_Exp2_vHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sig30 <- data_sig30[order(-data_sig30$significance),]
p3 <- ggplot(data_sig30,aes(Exp2_dHC_SW30_vs_Exp2_dHC_CON_ratio,Exp2_vHC_SW30_vs_Exp2_vHC_CON_ratio, color = -log10(data_sig30$significance))) + 
  geom_point() + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/max(-log10(data_sig30$significance)), 1)) +
  ggtitle("30 min Swim vs Control")+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(name = "dHC logFC Swim vs. Control") +
  scale_y_continuous(name = "vHC logFC Swim vs. Control") +
  labs(colour="-log10(FDR)") +
  theme_bw()

data_sig45 <- data_sig
data_sig45$significance <- apply(data_sig45[,c("Exp2_dHC_SW45_vs_Exp2_dHC_CON_FDR","Exp2_vHC_SW45_vs_Exp2_vHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sig45 <- data_sig45[order(-data_sig45$significance),]
p4 <- ggplot(data_sig45,aes(Exp2_dHC_SW45_vs_Exp2_dHC_CON_ratio,Exp2_vHC_SW45_vs_Exp2_vHC_CON_ratio, color = -log10(data_sig45$significance))) + 
  geom_point() + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/max(-log10(data_sig45$significance)), 1)) +
  ggtitle("45 min Swim vs Control")+
  scale_x_continuous(name = "dHC logFC Swim vs. Control") +
  scale_y_continuous(name = "vHC logFC Swim vs. Control") +
  labs(colour="-log10(FDR)") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()

plot_grid(p1,p2,p3,p4)


data_sigd6v15 <- data_sig
data_sigd6v15$significance <- apply(data_sigd6v15[,c("Exp1_dHC_SW06_vs_Exp1_dHC_CON_FDR","Exp2_dHC_SW15_vs_Exp2_dHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sigd6v15 <- data_sigd6v15[order(-data_sigd6v15$significance),]
pd6v15 <- ggplot(data_sigd6v15,aes(Exp1_dHC_SW06_vs_Exp1_dHC_CON_ratio,Exp2_dHC_SW15_vs_Exp2_dHC_CON_ratio, color = -log10(data_sigd6v15$significance))) + 
  geom_point() + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/max(-log10(data_sigd6v15$significance)), 1)) +
  ggtitle("Correlation swim effects 6 min vs 15 min in dHC")+
  scale_x_continuous(name = "dHC logFC Swim 6min vs. Control") +
  scale_y_continuous(name = "dHC logFC Swim 15min vs. Control") +
  labs(colour="-log10(FDR)") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()

data_sigd15v30 <- data_sig
data_sigd15v30$significance <- apply(data_sigd15v30[,c("Exp2_dHC_SW15_vs_Exp2_dHC_CON_FDR","Exp2_dHC_SW30_vs_Exp2_dHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sigd15v30 <- data_sigd15v30[order(-data_sigd15v30$significance),]
pd15v30 <- ggplot(data_sigd15v30,aes(Exp2_dHC_SW15_vs_Exp2_dHC_CON_ratio,Exp2_dHC_SW30_vs_Exp2_dHC_CON_ratio, color = -log10(data_sigd15v30$significance))) + 
  geom_point() + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/max(-log10(data_sigd15v30$significance)), 1)) +
  ggtitle("Correlation swim effects 15 min vs 30 min in dHC")+
  scale_x_continuous(name = "dHC logFC Swim 15min vs. Control") +
  scale_y_continuous(name = "dHC logFC Swim 30min vs. Control") +
  labs(colour="-log10(FDR)") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()

data_sigd30v45 <- data_sig
data_sigd30v45$significance <- apply(data_sigd30v45[,c("Exp2_dHC_SW30_vs_Exp2_dHC_CON_FDR","Exp2_dHC_SW45_vs_Exp2_dHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sigd30v45 <- data_sigd30v45[order(-data_sigd30v45$significance),]
pd30v45 <- ggplot(data_sigd30v45,aes(Exp2_dHC_SW30_vs_Exp2_dHC_CON_ratio,Exp2_dHC_SW45_vs_Exp2_dHC_CON_ratio, color = -log10(data_sigd30v45$significance))) + 
  geom_point() + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/max(-log10(data_sigd30v45$significance)), 1)) +
  ggtitle("Correlation swim effects 30 min vs 45 min in dHC")+
  scale_x_continuous(name = "dHC logFC Swim 30min vs. Control") +
  scale_y_continuous(name = "dHC logFC Swim 45min vs. Control") +
  labs(colour="-log10(FDR)") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()

data_sigv6v15 <- data_sig
data_sigv6v15$significance <- apply(data_sigv6v15[,c("Exp1_vHC_SW06_vs_Exp1_vHC_CON_FDR","Exp2_vHC_SW15_vs_Exp2_vHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sigv6v15 <- data_sigv6v15[order(-data_sigv6v15$significance),]
pv6v15 <- ggplot(data_sigv6v15,aes(Exp1_vHC_SW06_vs_Exp1_vHC_CON_ratio,Exp2_vHC_SW15_vs_Exp2_vHC_CON_ratio, color = -log10(data_sigv6v15$significance))) + 
  geom_point() + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/max(-log10(data_sigv6v15$significance)), 1)) +
  ggtitle("Correlation swim effects 6 min vs 15 min in vHC")+
  scale_x_continuous(name = "vHC logFC Swim 6min vs. Control") +
  scale_y_continuous(name = "vHC logFC Swim 15min vs. Control") +
  labs(colour="-log10(FDR)") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()

data_sigv15v30 <- data_sig
data_sigv15v30$significance <- apply(data_sigv15v30[,c("Exp2_vHC_SW15_vs_Exp2_vHC_CON_FDR","Exp2_vHC_SW30_vs_Exp2_vHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sigv15v30 <- data_sigv15v30[order(-data_sigv15v30$significance),]
pv15v30 <- ggplot(data_sigv15v30,aes(Exp2_vHC_SW15_vs_Exp2_vHC_CON_ratio,Exp2_vHC_SW30_vs_Exp2_vHC_CON_ratio, color = -log10(data_sigv15v30$significance))) + 
  geom_point() + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/max(-log10(data_sigv15v30$significance)), 1)) +
  ggtitle("Correlation swim effects 15 min vs 30 min in vHC")+
  scale_x_continuous(name = "vHC logFC Swim 15min vs. Control") +
  scale_y_continuous(name = "vHC logFC Swim 30min vs. Control") +
  labs(colour="-log10(FDR)") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()

data_sigv30v45 <- data_sig
data_sigv30v45$significance <- apply(data_sigv30v45[,c("Exp2_vHC_SW30_vs_Exp2_vHC_CON_FDR","Exp2_vHC_SW45_vs_Exp2_vHC_CON_FDR")],1,FUN = function(x){fisher(x)})
data_sigv30v45 <- data_sigv30v45[order(-data_sigv30v45$significance),]
pv30v45 <- ggplot(data_sigv30v45,aes(Exp2_vHC_SW30_vs_Exp2_vHC_CON_ratio,Exp2_vHC_SW45_vs_Exp2_vHC_CON_ratio, color = -log10(data_sigv30v45$significance))) + 
  geom_point() + 
  geom_density_2d() + 
  scale_colour_gradientn(colours = c("black", "red2", "gold2"),
                         values = c(0, -log10(0.05)/max(-log10(data_sigv30v45$significance)), 1)) +
  ggtitle("Correlation swim effects 30 min vs 45 min in vHC")+
  scale_x_continuous(name = "vHC logFC Swim 30min vs. Control") +
  scale_y_continuous(name = "vHC logFC Swim 45min vs. Control") +
  labs(colour="-log10(FDR)") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw()

plot_grid(pd6v15,pd15v30,pd30v45,pv6v15,pv15v30,pv30v45)

#Create Volcano plots 

data_volcano <- rbind(data.frame(logFC = data_results$Exp1_dHC_SW06_vs_Exp1_dHC_CON_ratio, PValue = data_results$Exp1_dHC_SW06_vs_Exp1_dHC_CON_p.val, FDR = data_results$Exp1_dHC_SW06_vs_Exp1_dHC_CON_FDR,Condition = "6min", Region = "dHC"),
                      data.frame(logFC = data_results$Exp1_vHC_SW06_vs_Exp1_vHC_CON_ratio, PValue = data_results$Exp1_vHC_SW06_vs_Exp1_vHC_CON_p.val, FDR = data_results$Exp1_vHC_SW06_vs_Exp1_vHC_CON_FDR,Condition = "6min", Region = "vHC"),
                      data.frame(logFC = data_results$Exp2_dHC_SW15_vs_Exp2_dHC_CON_ratio, PValue = data_results$Exp2_dHC_SW15_vs_Exp2_dHC_CON_p.val, FDR = data_results$Exp2_dHC_SW15_vs_Exp2_dHC_CON_FDR,Condition = "15min", Region = "dHC"),
                      data.frame(logFC = data_results$Exp2_vHC_SW15_vs_Exp2_vHC_CON_ratio, PValue = data_results$Exp2_vHC_SW15_vs_Exp2_vHC_CON_p.val, FDR = data_results$Exp2_vHC_SW15_vs_Exp2_vHC_CON_FDR,Condition = "15min", Region = "vHC"),
                      data.frame(logFC = data_results$Exp2_dHC_SW30_vs_Exp2_dHC_CON_ratio, PValue = data_results$Exp2_dHC_SW30_vs_Exp2_dHC_CON_p.val, FDR = data_results$Exp2_dHC_SW30_vs_Exp2_dHC_CON_FDR,Condition = "30min", Region = "dHC"),
                      data.frame(logFC = data_results$Exp2_vHC_SW30_vs_Exp2_vHC_CON_ratio, PValue = data_results$Exp2_vHC_SW30_vs_Exp2_vHC_CON_p.val, FDR = data_results$Exp2_vHC_SW30_vs_Exp2_vHC_CON_FDR,Condition = "30min", Region = "vHC"),
                      data.frame(logFC = data_results$Exp2_dHC_SW45_vs_Exp2_dHC_CON_ratio, PValue = data_results$Exp2_dHC_SW45_vs_Exp2_dHC_CON_p.val, FDR = data_results$Exp2_dHC_SW45_vs_Exp2_dHC_CON_FDR,Condition = "45min", Region = "dHC"),
                      data.frame(logFC = data_results$Exp2_vHC_SW45_vs_Exp2_vHC_CON_ratio, PValue = data_results$Exp2_vHC_SW45_vs_Exp2_vHC_CON_p.val, FDR = data_results$Exp2_vHC_SW45_vs_Exp2_vHC_CON_FDR,Condition = "45min", Region = "vHC"))


FCcutoff <- 0
fdrCutoff <- 0.05
fdrCutoff2 <- 0.5
change <- NULL
for(i in 1:dim(data_volcano)[1]){
  if(data_volcano[i,"logFC"] > FCcutoff & data_volcano[i,"FDR"] <= fdrCutoff){
    change<-append(change,"up-regulated 5% FDR")
  }
  else if(data_volcano[i,"logFC"] < -FCcutoff & data_volcano[i,"FDR"] <= fdrCutoff){
    change<-append(change,"down-regulated 5% FDR")
  }
  else if(data_volcano[i,"logFC"] > FCcutoff & data_volcano[i,"FDR"] <= fdrCutoff2){
    change<-append(change,"up-regulated 50% FDR")
  }
  else if(data_volcano[i,"logFC"] < -FCcutoff & data_volcano[i,"FDR"] <= fdrCutoff2){
    change<-append(change,"down-regulated 50% FDR")
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



ggplot(data_volcano,aes(logFC,PValue,colour=change))+
  geom_point()+
  scale_colour_manual(values = c("blue","#99ccff","grey30","red","salmon")) + 
  scale_y_continuous(trans=reverselog_trans(10)) + 
  theme_bw() +
  facet_grid(Region~Condition)


#Kmeans clustering

clustering_dat <- data_sig[,c("Exp1_dHC_SW06_vs_Exp1_dHC_CON_ratio","Exp1_vHC_SW06_vs_Exp1_vHC_CON_ratio","Exp2_dHC_SW15_vs_Exp2_dHC_CON_ratio","Exp2_vHC_SW15_vs_Exp2_vHC_CON_ratio","Exp2_dHC_SW30_vs_Exp2_dHC_CON_ratio","Exp2_vHC_SW30_vs_Exp2_vHC_CON_ratio","Exp2_dHC_SW45_vs_Exp2_dHC_CON_ratio","Exp2_vHC_SW45_vs_Exp2_vHC_CON_ratio")]
row.names(clustering_dat) <- data_sig$name

Results <- kmeans(clustering_dat,centers = 5, nstart = 20)
clustering_dat <- data.frame(clustering_dat, cluster = Results$cluster)
clustering_dat$name <- row.names(clustering_dat)

clustering_dat <- gather(clustering_dat,key = "Comparison", value = log2FC, -cluster, -name)
clustering_dat$TimePointNumeric <-  sapply(clustering_dat$Comparison, FUN = function(x){switch(x,
                                                                                               Exp1_dHC_SW06_vs_Exp1_dHC_CON_ratio = as.numeric(6),
                                                                                               Exp1_vHC_SW06_vs_Exp1_vHC_CON_ratio = as.numeric(6),
                                                                                               Exp2_dHC_SW15_vs_Exp2_dHC_CON_ratio = as.numeric(15),
                                                                                               Exp2_vHC_SW15_vs_Exp2_vHC_CON_ratio = as.numeric(15),
                                                                                               Exp2_dHC_SW30_vs_Exp2_dHC_CON_ratio = as.numeric(30),
                                                                                               Exp2_vHC_SW30_vs_Exp2_vHC_CON_ratio = as.numeric(30),
                                                                                               Exp2_dHC_SW45_vs_Exp2_dHC_CON_ratio = as.numeric(45),
                                                                                               Exp2_vHC_SW45_vs_Exp2_vHC_CON_ratio = as.numeric(45))})
clustering_dat$Region <-  sapply(clustering_dat$Comparison, FUN = function(x){switch(x,
                                                                                               Exp1_dHC_SW06_vs_Exp1_dHC_CON_ratio = "dHC",
                                                                                               Exp1_vHC_SW06_vs_Exp1_vHC_CON_ratio = "vHC",
                                                                                               Exp2_dHC_SW15_vs_Exp2_dHC_CON_ratio = "dHC",
                                                                                               Exp2_vHC_SW15_vs_Exp2_vHC_CON_ratio = "vHC",
                                                                                               Exp2_dHC_SW30_vs_Exp2_dHC_CON_ratio = "dHC",
                                                                                               Exp2_vHC_SW30_vs_Exp2_vHC_CON_ratio = "vHC",
                                                                                               Exp2_dHC_SW45_vs_Exp2_dHC_CON_ratio = "dHC",
                                                                                               Exp2_vHC_SW45_vs_Exp2_vHC_CON_ratio = "vHC")})

ggplot(clustering_dat,aes(TimePointNumeric,log2FC, group = name, color = as.factor(cluster))) + geom_line() + facet_wrap(vars(cluster, Region), scales = "free_y", ncol = 10) + theme_bw()
assign <- match(clustering_dat$name,data_results$name)
clustering_dat$Protein <- data_results[assign,"Protein"]
clustering_dat$GeneSymbol <- data_results[assign,"GeneSymbol"]

#write.table(clustering_dat,"Clusterdat_n25.csv", sep = ";", row.names = F)
#write.table(clustering_dat,"Clusterdat_n10.csv", sep = ";", row.names = F)
#write.table(clustering_dat,"Clusterdat_n5.csv", sep = ";", row.names = F)



#LIneplots for subset of targets
targets <- c("_RAPS[Phospho_STY_NL_98 (STY)]PVVS[Phospho_STY_NL_98 (STY)]PTELSK_")

plotdat <- data_volcano[data_volcano$name %in% targets,]
plotdat$name <- factor(plotdat$name, levels = targets)
ggplot(plotdat,aes(TimePointNumeric,logFC, color = Region)) + geom_path() + scale_x_continuous(breaks = c(6,15,30,45)) + geom_hline(yintercept = 0) + facet_grid(name~.) + theme_bw()




data_volcano$TimePointNumeric <-  sapply(data_volcano$Condition, FUN = function(x){switch(x,
                                                                                               '6min' = as.numeric(6),
                                                                                               '15min' = as.numeric(15),
                                                                                               '30min' = as.numeric(30),
                                                                                               '45min' = as.numeric(45))})



targets <- c("_RAPS[Phospho_STY_NL_98 (STY)]PVVS[Phospho_STY_NL_98 (STY)]PTELSK_")

plotdat <- matr["_RAPS[Phospho_STY_NL_98 (STY)]PVVS[Phospho_STY_NL_98 (STY)]PTELSK_",]
assign <- match(names(plotdat),s2c$FileName)
plotdat <- data.frame(Int = plotdat, Hemisphere = s2c[assign,"Region"], condition = s2c[assign,"ShortCond"])

#dHC vs vHC at 6min control
data <- read.table("20191031_174447_20191028_AllPhosphoData_Report.csv",sep = ",", header = T)
s2c <- read.table("metadata/s2c_allphos.csv", sep = ";", header = T)
rownames(s2c) <- s2c$FileName

length(unique(data$EG.PrecursorId))

matr <- data[,c("R.FileName","EG.ModifiedPeptide","FG.Quantity")]
matr <- aggregate(FG.Quantity ~ R.FileName + EG.ModifiedPeptide, data = matr, sum)


matr <- spread(matr, key = R.FileName, value = FG.Quantity)

rownames(matr) <- matr[,1]
matr <- matr[,-1]


s2csub <- s2c
mat <- matr[,names(matr) %in% s2csub$FileName]
mat$Protein.names <- row.names(mat)
mat$Protein.IDs <- paste("ID_",c(1:nrow(mat)),sep = "")
data_unique <- make_unique(mat, "Protein.names", "Protein.IDs", delim = ";")

experimental_design <- data.frame(label = as.character(s2csub[,"FileName"]), 
                                  condition = as.character(paste(s2csub[,"ShortExp"],s2csub[,"Region"],s2csub[,"ShortCond"], sep = "_")), 
                                  replicate = as.integer(s2csub[,"Replicate"]))
experimental_design[,c("label","condition")] <- lapply(experimental_design[,c("label","condition")], as.character)


#Create DEP data
LFQ_columns <- na.omit(match(s2csub$FileName,colnames(data_unique)))
data_se <- make_se(data_unique, LFQ_columns, experimental_design)
data_filt <- filter_missval(data_se, thr = 2)
data_norm <- normalize_vsn(data_filt)
assays(data_filt)$vsn <- assay(data_norm)
data_imp <- impute(data_filt, fun = "MinProb", q = 0.01)
assays(data_filt)$imputed <- assay(data_imp)

data_diff_manual <- test_diff(data_imp, type = "manual", test = c("Exp1_dHC_CON_vs_Exp1_vHC_CON"))

#Export and build Results data.frame
data_results <- get_results(dep_man)
assign <- match(data_results$name,data$EG.ModifiedPeptide)
data_results$Protein <- data[assign,"PG.ProteinAccessions"]
data_results$ModifiedSequence <- data[assign,"EG.PTMLocalizationProbabilities"]
data_results$StrippedSequence <- data[assign,"EG.StrippedSequence"]
data_results$Positions <- data[assign,"EG.PTMPositions..Phospho_STY_NL_98..STY.."]
data_results$Probabilities <- data[assign,"EG.PTMProbabilities..Phospho_STY_NL_98..STY.."]

data_results$Exp1_dHC_CON_vs_Exp1_vHC_CON_FDR <- p.adjust(data_results$Exp1_dHC_CON_vs_Exp1_vHC_CON_p.val, method = "fdr")
data_volcano <- data.frame(FDR = data_results$Exp1_dHC_CON_vs_Exp1_vHC_CON_FDR, PValue = data_results$Exp1_dHC_CON_vs_Exp1_vHC_CON_p.val, logFC = data_results$Exp1_dHC_CON_vs_Exp1_vHC_CON_ratio)


FCcutoff <- 0
fdrCutoff <- 0.05
change <- NULL
for(i in 1:dim(data_volcano)[1]){
  if(data_volcano[i,"logFC"] > FCcutoff & data_volcano[i,"FDR"] <= fdrCutoff){
    change<-append(change,"up-regulated 5% FDR")
  }
  else if(data_volcano[i,"logFC"] < -FCcutoff & data_volcano[i,"FDR"] <= fdrCutoff){
    change<-append(change,"down-regulated 5% FDR")
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



ggplot(data_volcano,aes(logFC,PValue,colour=change))+
  geom_point()+
  scale_colour_manual(values = c("blue","black","red")) + 
  scale_y_continuous(trans=reverselog_trans(10)) + 
  theme_bw()
