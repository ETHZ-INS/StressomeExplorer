library(shiny)
library(SummarizedExperiment)
library(SEtools)
library(ggplot2)
library(cowplot)

## Toggle this to enable/disable password protection
Logged = TRUE;

PASSWORD <- data.frame(username=c("gbohacek"),password=c("letmein"))

ff <- list.files("../data/seq/", pattern="\\.SE\\.", full.names=TRUE)
names(ff) <- gsub("\\.SE\\.rds","",basename(ff))
SEs <- lapply(ff,FUN=readRDS)
phos <- readRDS("../data/phos/PhosphoData.SE.rds")
prot <- readRDS("../data/prot/4HFST_LFQ.SE.rds")

tgl <- list()

formatDEA <- function(x){
  for(f in colnames(x)){
    if(is.numeric(x[[f]])) x[[f]] <- plgINS::dround(x[[f]])
  }
  x
}

grepGene <- function(x,g){
  g <- paste0("^",g,"\\.|^",g,"$")
  g <- unlist(sapply(g,FUN=function(i) grep(i, row.names(x))))
  x[g,,drop=FALSE]
}

shinyServer(function(input, output, session) {
  source("Login.R", local = TRUE)
  g <- sapply(strsplit(unlist(lapply(SEs, row.names)),".",fixed=T),FUN=function(x) x[1])
  updateSelectizeInput(session, "gene_input", choices=sort(unique(g)),selected = "Fos", server=T)

  observe({
    if (USER$Logged == TRUE) {
      output$sidebarui <- renderUI(list(
        sidebarMenu(
          menuItem("Gene", tabName="tab_gene")
        )))
      }
    updateSelectInput(session, "select_datasets", choices=names(SEs))
  })


  PlotHeight <- reactive(
    return(100 +100 * sum(rowData(phos)$GeneSymbol == input$gene_input,na.rm = T))
  )

  ############
  ### START GENE PLOT
  output$gene_name <- renderText({input$gene_input})

  output$gene_plot <- renderPlot({
    ll <- lapply(SEs, FUN=function(x){
      g <- row.names(grepGene(x,input$gene_input))
      x <- SEtools::meltSE(x,g,assayName=intersect(c("log2FC","logFC","corrected","logcpm","lognorm"), assayNames(x))[1])
      colnames(x)[ncol(x)] <- "value"
      x
    })
    ll <- ll[sapply(ll,nrow)>0]
    d <- as.data.frame(data.table::rbindlist(ll, fill=TRUE, idcol="Dataset"))
    d <- d[,!duplicated(colnames(d))]
    
    if(!input$select_logaxis){
      d$value <- exp(d$value)
    }
    
    p <- list()
    TS <- d[d$Experiment == "TimeSeriesFST",]
    if(nrow(TS) > 0){
    p[["TS"]] <- ggplot(TS, aes(Condition2, value, color = Condition2))
    p[["TS"]] <- p[["TS"]] + 
      facet_grid(.~Region) + 
      ggtitle(paste(input$gene_input, " in the acute stress time series", sep = "")) +
      scale_color_manual(values = c("#441C53","#43488A","#277A91","#21A885","#81C456","#F7E820"))
    }else{p[["TS"]] <- ggplot() + ggtitle(paste(input$gene_input, " not present in the acute stress time series", sep = ""))}
    
    sex <- d[d$Experiment == "FSTFemale45min4h",]
    if(nrow(sex) > 0){
    p[["sex"]] <- ggplot(sex, aes(Condition2, value, color = Condition2))
    p[["sex"]] <- p[["sex"]] + 
      facet_grid(Region~Sex) + 
      ggtitle(paste(input$gene_input, " in males and females after AS (vHC)", sep = "")) + 
      scale_color_manual(values = c("#441C53","#43488A","#F7E820"))
    }else{p[["sex"]] <- ggplot() + ggtitle(paste(input$gene_input, " not present in males and females after AS (vHC)", sep = ""))}
    
    hem <- d[d$Experiment == "LeftvsRight",]
    if(nrow(hem) > 0){
    p[["hem"]] <- ggplot(hem, aes(Condition2, value, color = Condition2))
    p[["hem"]] <- p[["hem"]] + 
      facet_grid(Hemisphere~Region) + 
      ggtitle(paste(input$gene_input, " in left and right hemispheres after AS", sep = "")) +
      scale_color_manual(values = c("#441C53","#43488A"))
    }else{p[["hem"]] <- ggplot() + ggtitle(paste(input$gene_input, " not present in left and right hemispheres after AS", sep = ""))}
    
    
    CMV <- d[d$Experiment == "CMV-TRAP",]
    if(nrow(CMV) > 0){
    CMV$protocol <- ifelse(CMV$SampleProcessing == "WholeRNA","Pre-IP","Post-IP")
    p[["CMV"]] <- ggplot(CMV, aes(Condition2, value, color = Condition2))
    p[["CMV"]] <- p[["CMV"]] + 
      facet_grid(Region~protocol) + 
      ggtitle(paste(input$gene_input, " in CMV-nuTRAP (ubiquitous expression) after AS (vHC)", sep = ""))+
      scale_color_manual(values = c("#441C53","#43488A"))
    }else{p[["CMV"]] <- ggplot() + ggtitle(paste(input$gene_input, " not present in CMV-nuTRAP (ubiquitous expression) after AS (vHC)", sep = ""))}
    
    
    CAMK2A <- d[d$Experiment %in% c("CAMK2A-TRAP-2-2M","CAMK2A-TRAP","VIAAT-TRAP"),]
    if(nrow(CAMK2A) > 0){
    CAMK2A$Type <- ifelse(CAMK2A$Experiment == "VIAAT-TRAP","VIAAT-bacTRAP","CAMK2A-nuTRAP")
    p[["CAMK2A"]] <- ggplot(CAMK2A, aes(Condition2, value, color = Condition2))
    p[["CAMK2A"]] <- p[["CAMK2A"]] + 
      facet_grid(Type~Region, scales = "free_x") + 
      ggtitle(paste(input$gene_input, " in excitatory (CAMK2A) and inhibitory (VIAAT) neurons after AS", sep = "")) +
      scale_color_manual(values = c("#441C53","#43488A","#277A91"))
    }else{p[["CAMK2A"]] <- ggplot() + ggtitle(paste(input$gene_input, " not present in excitatory (CAMK2A) and inhibitory (VIAAT) neurons after AS", sep = ""))}
    
    
    
    for(i in names(p)){
      yaxname<- ifelse(input$select_logaxis,"logcpm","cpm")
      p[[i]] <- p[[i]] + xlab("") + ylab(yaxname)
    if(input$select_plottype == "violin plot"){
      p[[i]] <- p[[i]] + geom_violin()
    }else if(input$select_plottype == "box plot"){
      p[[i]] <- p[[i]] + geom_boxplot(outlier.shape = NA)
    }
    if(input$select_plotpoints){
      p[[i]] <-p[[i]] + geom_point(position = position_jitterdodge())
    }
    p[[i]] <- p[[i]] + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank())
    }
    
    print(plot_grid(p[["TS"]], p[["sex"]], p[["hem"]],p[["CMV"]],p[["CAMK2A"]],ncol = 1, rel_heights = c(1,1,1.5,1,1.5)))
  })
  
  
  output$gene_inList <- renderPrint({
    if(is.null(input$gene_input) || length(tgl)==0) return(NULL)
    x <- which(sapply(tgl,g=input$gene_input,FUN=function(x,g){ g %in% x}))
    if(length(x)==0) return("This gene is included in none of the registered genelists.")
    cat("This gene is included in the following list(s):\n\n")
    for(i in names(tgl)[x]) cat(paste(i,"\n"))
  })
  ### END GENE PLOT
  ############

  
  ### PHOSPHOS PLOT
  ############
  output$phospho_plot <- renderPlot(height = function(){PlotHeight()},{

    subset <- phos[which(rowData(phos)$GeneSymbol == input$gene_input),]
    
    if(nrow(subset) > 0){
    subset <- subset[order(rowData(subset)$EarliestPos),]
    rowData(subset)$pepNr <- paste("Pep_", 1:nrow(subset),sep = "")
    rowData(subset)$pepNr <- factor(rowData(subset)$pepNr, levels = rowData(subset)$pepNr)
    rdat <- rowData(subset)
    pepdats <- NULL
    for(i in rownames(rdat)){
      sq <- as.character(rdat[i,"StrippedSequence"])
      po <- rdat[i,"EarliestPos"]:(rdat[i,"EarliestPos"]+nchar(sq) - 1)
      prob <- rep(0,nchar(sq))
      if(rdat[i,"Probabilities"] != ""){
        anot <- as.numeric(unlist(strsplit(as.character(rdat[i,"Positions"]), split = ";")))
        anotprob <- as.numeric(unlist(strsplit(as.character(rdat[i,"Probabilities"]), split = ";")))
        prob[anot] <- anotprob
      }
      pepdats <- rbind(pepdats, data.frame(ID = rdat[i,"pepNr"], name = i, AA = unlist(strsplit(sq, split = "")), Position = po, Probability = prob))
    }
    
    
    p0 <- ggplot(pepdats,aes(Position, "AA", color = "detected sequence")) + 
      scale_color_manual(values = "darkolivegreen") + 
      geom_point() + 
      ggtitle("Assay coverage") + 
      theme_bw()
    p1 <- ggplot(pepdats,aes(Position, "AA", label = AA, color = Probability)) + 
      geom_text() + 
      facet_wrap(~ID, scales = "free", ncol = 1) + 
      theme_bw() + 
      ylab("") + 
      scale_color_gradient(low = "black", high = "steelblue") + 
      ggtitle("Phosphorylation Probabilities") +
      theme(legend.position="top")
    
    cdat <- SEtools::meltSE(subset,names(subset),assayName=intersect(c("vsn","imputed"), assayNames(subset))[1])
    assign <- match(cdat$feature,rownames(rowData(subset)))
    cdat$ID <- rowData(subset)$pepNr[assign]
    cdat$TimePoint <- factor(cdat$TimePoint, levels = c("0min","6min","15min","30min","45min"))
    
    if(input$select_logaxis){
      yaxname <- "log2 Intensity"
    }else{
      cdat$vsn <- 2^cdat$vsn
      yaxname <- "Intensity"
    }
    
    p2 <- ggplot(cdat[cdat$Region == "vHC",],aes(TimePoint,vsn, color = TimePoint)) + 
      facet_grid(ID~Experiment, scales = "free") + 
      theme_bw() + ggtitle("vHC") + 
      theme(legend.position = "top", legend.title = element_blank() ,axis.text.x = element_text(angle = 45, hjust = 1)) + 
      scale_color_manual(values = c("#441C53","#24798F","#25A885","#81C456","#F8E716")) +
      ylab(yaxname)
    
    p3 <- ggplot(cdat[cdat$Region == "dHC",],aes(TimePoint,vsn, color = TimePoint)) +
      facet_grid(ID~Experiment, scales = "free") + 
      theme_bw() + ggtitle("dHC") + 
      theme(legend.position = "top", legend.title = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1)) + 
      scale_color_manual(values = c("#441C53","#24798F","#25A885","#81C456","#F8E716")) +
      ylab(yaxname)
    
    
    if(input$select_plottype == "violin plot"){
      p2 <- p2 + geom_violin()
      p3 <- p3 + geom_violin()
    }else if(input$select_plottype == "box plot"){
      p2 <- p2 + geom_boxplot(outlier.shape = NA)
      p3 <- p3 + geom_boxplot(outlier.shape = NA)
    }
    if(input$select_plotpoints){
      p2 <-p2 + geom_point(position = position_jitterdodge())
      p3 <-p3 + geom_point(position = position_jitterdodge())
    }
    
    p4 <- plot_grid(p1,p2,p3, nrow = 1, rel_widths = c(1.7,1,1))
    plot_grid(p0,p4,ncol = 1, rel_heights = c(1,nrow(subset)))
    }
    else
    {
      p1 <- ggplot() + ggtitle(paste("No phosphorylation data for ", input$gene_input, sep = "")) + theme_bw()
      print(p1)
    }
  })
  ### END PHOSPHO PLOT
  ############
  
  
  ### LFQ PLOT
  ############
  output$protplot_plot <- renderPlot({
    
    subset <- prot[which(rowData(prot)$Genes == input$gene_input),]
    
    if(nrow(subset) == 1){
    cdat <- SEtools::meltSE(subset,names(subset),assayName=intersect(c("lognorm.imputed","lognorm","intensity"), assayNames(subset))[2])
    cdat$Condition <- ifelse(cdat$Condition == "Swim","Swim 4h", "Homecage")
    
    if(input$select_logaxis){
      yaxname <- "log2 Intensity"
    }else{
      cdat$lognorm <- 2^cdat$lognorm
      yaxname <- "Intensity"
    }
    
    p2 <- ggplot(cdat,aes(1,lognorm, color = SubRegion)) + 
      facet_grid(.~Region) + 
      theme_bw() + 
      xlab("") +
      scale_color_manual(values = c("#CEB733","#34489C","#12A04D")) +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
      ggtitle(paste(input$gene_input, " in regions and sub-regions (pooled groups)", sep = "")) + ylab(yaxname)
    
    p3 <- ggplot(cdat,aes(Condition,lognorm, color = Condition)) + 
      scale_color_manual(values = c("#441C53","#F8E716")) +
      facet_wrap(Region~SubRegion, scales = "free") + 
      theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(paste(input$gene_input, " Controls vs. 4h AS", sep = "")) + ylab(yaxname)
    
    if(input$select_plottype == "violin plot"){
      p2 <- p2 + geom_violin()
      p3 <- p3 + geom_violin()
    }else if(input$select_plottype == "box plot"){
      p2 <- p2 + geom_boxplot(outlier.shape = NA)
      p3 <- p3 + geom_boxplot(outlier.shape = NA)
    }
    if(input$select_plotpoints){
      p2 <-p2 + geom_point(position = position_jitterdodge())
      p3 <-p3 + geom_point(position = position_jitterdodge())
    }
    plot_grid(p2,p3, ncol = 1)
    }else{
      p1 <- ggplot() + ggtitle(paste("No proteomic data for ", input$gene_input, sep = ""))
      print(p1)
    }
  })
  ### END LFQ PLOT
  ############
})

