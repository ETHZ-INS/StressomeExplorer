library(shiny)
library(SummarizedExperiment)
library(SEtools)
library(ggplot2)
library(cowplot)

data_dir <- "../data"
fig_dir <- "../fig"

ff <- list.files(file.path(data_dir,"seq"), pattern="\\.SE\\.", full.names=TRUE)
names(ff) <- gsub("\\.SE\\.rds","",basename(ff))
SEs <- lapply(ff,FUN=readRDS)
phos <- readRDS(file.path(data_dir,"phos","PhosphoData.SE.rds"))
prot <- readRDS(file.path(data_dir, "prot", "4HFST_LFQ.SE.rds"))
SN <- readRDS(file.path(data_dir, "snRNAseq", "snRNAseq.SE.rds"))

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

phosphoHeat <- function(subset, position_width=3, assayName="vsn", ...){
  library(gridtext)
  library(ComplexHeatmap)
  subset <- subset[!grepl("Ambiguous",rowData(subset)$PhosphoSites),]
  subset <- subset[order(rowData(subset)$EarliestPos),]
  rdat <- rowData(subset)
  sqs <- strsplit(as.character(rdat$StrippedSequence),"")
  pos <- lapply(strsplit(as.character(rdat$Positions),";"), as.integer)
  probs <- lapply(strsplit(as.character(rdat$Probabilities),";"), FUN=function(x){
    as.integer(as.numeric(x)*100)
  })
  cols <- viridisLite::cividis(101)
  sq <- mapply(s=sqs, po=pos, prob=probs, FUN=function( s, po, prob){
    if(length(w <- which(prob>0))>0){
      s[po[w]] <- paste0("<span style='color: ", cols[1L+prob[w]],"'>", s[po[w]], "</span>")
    }
    paste(s,collapse="")
  })
  m <- matrix(0L, ncol=max(rdat$EarliestPos+lengths(sqs)), nrow=length(sq))
  for(i in seq_len(nrow(rdat))){
    m[i,rdat$EarliestPos[i]:(rdat$EarliestPos[i]+length(sqs[[i]])-1)] <- 1L
  }
  h1 <- Heatmap(m, cluster_rows=FALSE, cluster_columns=FALSE,
                col=c("lightgrey","darkblue"),
                column_title="Position", show_heatmap_legend=FALSE,
                heatmap_width=unit(position_width, "cm"))
  ra <- rowAnnotation(text=anno_text(gt_render(sq), gp=gpar(fontfamily="mono", fontface="bold")))
  h2 <- sechm(subset[,order(subset$Region,subset$TimePoint)], genes=row.names(subset),
              do.scale=FALSE, gaps_at=c("Region"), anno_columns=c("TimePoint","Region"),
              show_rownames=FALSE, sortRowsOn=NULL, cluster_rows=FALSE, assayName=assayName,
              right_annotation=ra, breaks=0.99, ...)
  draw(h1+h2, annotation_legend_list = list(
    Legend(col_fun=circlize::colorRamp2(breaks=c(0:100), colors=cols),
           title = "Phosphorylation\nprobability", at=c(0,50,100))), merge_legends=TRUE)
}


shinyServer(function(input, output, session) {
  g <- sapply(strsplit(unlist(lapply(SEs, row.names)),".",fixed=T),FUN=function(x) x[1])
  updateSelectizeInput(session, "gene_input", choices=sort(unique(g)),selected = "Fos", server=T)

  PlotHeight_Phos <- reactive(
    return(100 + 10 * sum(rowData(phos)$GeneSymbol == input$gene_input,na.rm = T))
  )

  PlotHeight_Phospep <- reactive(
    return(100 +100 * sum(rowData(phos)$GeneSymbol == input$gene_input,na.rm = T))
  )
  
  ############
  ### START Experimental design images
output$EDTS <- renderImage({
    list(src = paste0(fig_dir,"/Figure3A_EXPERIMENTALDESIGN.png"),width = 290,
         height = 190)
  }, deleteFile = FALSE)
  output$FvM <- renderImage({
    list(src = paste0(fig_dir,"/FigureFvM_EXPERIMENTALDESIGN.png"),width = 380,
         height = 100)
  }, deleteFile = FALSE)
  output$LvR <- renderImage({
    list(src = paste0(fig_dir,"/FigureLvR_EXPERIMENTALDESIGN.png"),width = 400,
         height = 80)
  }, deleteFile = FALSE)
  output$CMV <- renderImage({
    list(src = paste0(fig_dir,"/FigureCMV_EXPERIMENTALDESIGN.png"),width = 170,
         height = 80)
  }, deleteFile = FALSE)
  output$TRAP <- renderImage({
    list(src = paste0(fig_dir,"/FigureTRAPs_EXPERIMENTALDESIGN.png"),width = 400,
         height = 200)
  }, deleteFile = FALSE)
  output$EDphos <- renderImage({
    list(src = paste0(fig_dir,"/Figure2A_EXPERIMENTALDESIGN.png"),width = 360,
         height = 200)
  }, deleteFile = FALSE)
  output$EDphos2 <- renderImage({
    list(src = paste0(fig_dir,"/Figure2A_EXPERIMENTALDESIGN.png"),width = 360,
         height = 200)
  }, deleteFile = FALSE)
  output$EDprot <- renderImage({
    list(src = paste0(fig_dir,"/Figure7A_EXPERIMENTALDESIGN.png"),width = 400,
         height = 90)
  }, deleteFile = FALSE)
  
  ### END Experimental design images
  ############ 

  ############
  ### START GENE PLOT
  output$gene_name <- renderText({input$gene_input})

  output$availability <- renderText({
    snRNAdat <- SEtools::meltSE(SN, input$gene_input, assayName="logcpm")
    s1 <- paste("<b>availablility of", input$gene_input,":</b>" , sep = " ")
    s2 <- "<font color=\"#00CC00\">Transcriptome data</font>"
    s3 <- ifelse(nrow(snRNAdat) > 0, "<font color=\"#00CC00\">snRNA data</font>", "<font color=\"#FF0000\">no snRNA data</font>")
    s4 <- ifelse(sum(rowData(prot)$Genes == input$gene_input, na.rm = T) > 0, "<font color=\"#00CC00\">Proteome data</font>", "<font color=\"#FF0000\">no Proteome data</font>")
    s5 <- ifelse(sum(rowData(phos)$GeneSymbol == input$gene_input, na.rm = T) > 0, "<font color=\"#00CC00\">Phospho data</font>", "<font color=\"#FF0000\">no Phospho data</font>")
    paste(s1,s2,s3,s4,s5, sep = "<br/>")
  })

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
  output$phospho_plot <- renderPlot(height = function(){PlotHeight_Phos()},{

    subset <- phos[which(rowData(phos)$GeneSymbol == input$gene_input),]

    if(nrow(subset) > 0){
      subset <- subset[order(rowData(subset)$EarliestPos),]
      subset$TimePoint <- factor(subset$TimePoint, levels = c("0min","6min","15min","30min","45min"))
      ac <- list(TimePoint = c("0min"="#441C53","6min"="#24798F","15min"="#25A885","30min"="#81C456","45min"="#F8E716"),
                 Region = c("dHC"="#B29249","vHC"="#8A9DAF"))
      phosphoHeat(subset, anno_colors = ac, assayName=input$phos_assay)
    }
    else
    {
      p1 <- ggplot() + ggtitle(paste("No phosphorylation data for ", input$gene_input, sep = "")) + theme_bw()
      print(p1)
    }
  })

  output$phospho_pep_plot <- renderPlot(height = function(){PlotHeight_Phospep()},{

    subset <- phos[which(rowData(phos)$GeneSymbol == input$gene_input),]
    
    int_breaks <- function(x, n = 10) {
      l <- pretty(x, n)
      l[abs(l %% 1) < .Machine$double.eps ^ 0.5] 
    }

    if(nrow(subset) > 0){
    subset <- subset[!grepl("Ambiguous",rowData(subset)$PhosphoSites),]
    subset <- subset[order(rowData(subset)$EarliestPos),]
    rowData(subset)$pepNr <- paste("Pep_", 1:nrow(subset),sep = "")
    rowData(subset)$pepNr <- factor(rowData(subset)$pepNr, levels = rowData(subset)$pepNr)
    rdat <- rowData(subset)

    sq <- as.character(rdat$StrippedSequence)
    pos <- lapply(strsplit(as.character(rdat$Positions),";"), as.integer)
    probs <- mapply(ep=rdat$EarliestPos, nc=nchar(sq), pos=pos,
           prob=strsplit(as.character(rdat$Probabilities), ";"),
           FUN=function(ep,nc,pos,prob){
             p <- rep(0,nc)
             p[as.integer(pos)] <- as.numeric(prob)
             p
           })
    pepdats <- data.frame(
      ID=rep(rdat$pepNr, lengths(probs)),
      name=rep(row.names(rdat), lengths(probs)),
      AA=unlist(strsplit(sq,"")),
      Position=unlist(mapply(a=as.integer(rdat$EarliestPos),
                             b=lengths(probs), FUN=function(a,b) a:(a+b-1))),
      Probability=unlist(probs)
    )

    p0 <- ggplot(pepdats,aes(Position, "AA", color = "detected sequence")) +
      scale_color_manual(values = "darkolivegreen") +
      geom_point() +
      ggtitle("Assay coverage") +
      theme_bw()
    p1 <- ggplot(pepdats,aes(Position, "AA", label = AA, color = Probability)) +
      geom_text(fontface='bold') +
      facet_wrap(~ID, scales = "free", ncol = 1) +
      theme_bw() +
      ylab("") +
      scale_color_gradient(low = "#000033", high = "#FF9900") +
      scale_x_continuous(breaks = int_breaks) +
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

  output$snrna_plot <- renderPlot({
    cdat <- SEtools::meltSE(SN, input$gene_input, assayName="logcpm")
    if(nrow(cdat) > 0){
    levels(cdat$condition) <- c("Homecage", "Swim 45min")
    p1 <- ggplot(cdat, aes(cluster, logcpm, color = cluster)) + geom_point(position = position_jitterdodge()) +
      theme_bw() +
      theme(legend.position="none") +
      labs(x="", y="sctransform-normalized expression") +
      ggtitle("Expression in cell types (pooled groups)")
    p2 <- ggplot(cdat, aes(condition, logcpm, color=condition)) + geom_point(position = position_jitterdodge()) + facet_wrap(~cluster, scales="free_y") +
      scale_color_manual(values = c("#441C53","#24798F")) +
      theme_bw() +
      theme(legend.position="none") +
      labs(x="", y="sctransform-normalized expression") +
      ggtitle("Expression after acute stress")
    
    if(input$select_plottype == "violin plot"){
      p1 <- p1 + geom_violin()
    }else if(input$select_plottype == "box plot"){
      p1 <- p1+ geom_boxplot(outlier.shape = NA)
    }
    if(input$select_plotpoints){
      p1 <-p1 + geom_point(position = position_jitterdodge())
    }
    
    plot_grid(p1,p2,ncol = 1, rel_heights = c(1,3))
    }else{
      ggplot() + ggtitle(paste("No snRNA data for ", input$gene_input, sep = ""))
    }

  })

})

