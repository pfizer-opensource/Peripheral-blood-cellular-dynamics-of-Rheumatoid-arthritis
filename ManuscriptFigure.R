#Load Libraries
library(ggpubr)
library(ggplot2)
library(cowplot)

#Load Files
  gset     <- readRDS("~/data_objects/CR uploads_KI_delivery_combine_delivery_data_objects_eset_qced.rds")
  Cellgset <- readRDS("~/data_objects/CR uploads_KI_delivery_combine_delivery_data_objects_contributions_qced.rds")
  ctDiffMTX <- read.table("~/ct_diff/CR uploads_KI_delivery_combine_delivery_ct_diff_A_ct_diff_delivery.csv",
                       sep=",", header=T)
  ctDiffTNF <- read.table("~/ct_diff/CR uploads_KI_delivery_combine_delivery_ct_diff_A_ct_diff_delivery.csv",
                       sep=",", header=T)
  gtDiffMTX <- read.table("~/combine_delivery/gx_diff/CR uploads_KI_delivery_combine_delivery_gx_diff_A_all_results_delivery.csv",
                          sep=",", header=T)
  gtDiffTNF <- read.table("~/sjelinsky/combine_delivery/gx_diff/CR uploads_KI_delivery_combine_delivery_gx_diff_B_all_results_delivery.csv",
                          sep=",", header=T)


    ctDiffMTX$Treatment <- "MTX"
    ctDiffTNF$Treatment <- "TNF"
    
    ctDiffCombine <- rbind(ctDiffMTX, ctDiffTNF)

    cellorder <- c("monocyte", "CD14-positive, CD16-positive monocyte", "CD14-positive, CD16-negative classical monocyte", "CD14-low, CD16-positive monocyte",                                      
               "myeloid dendritic cell", "dendritic cell" ,  "plasmacytoid dendritic cell" ,
                  "T-helper 17 cell" ,  "T-helper 1 cell", "central memory CD8-positive, alpha-beta T cell", "CD4-positive, alpha-beta T cell",                            
                  "CD8-positive, alpha-beta T cell",           
                  "effector memory CD4-positive, alpha-beta T cell",            "effector memory CD8-positive, alpha-beta T cell",           
                  "effector memory RA CD8-positive, alpha-beta T cell (TEMRA)" ,                                           
                  "mature NK T cell" ,                     
                  "naive thymus-derived CD4-positive, alpha-beta T cell",       "naive thymus-derived CD8-positive, alpha-beta T cell",      
                  "regulatory T cell" ,                    
                  "T-helper 2 cell",
                  "central memory CD4-positive, alpha-beta T cell", 
                  "CD16-positive, CD56-dim natural killer cell", "natural killer cell", "mature natural killer cell",
                  "plasma cell", "memory B cell" , "class switched memory B cell", "mature B cell" , "naive B cell", 
                  "eosinophil", "granulocyte","neutrophil")
###### 
#Figure 1
######    
#Figure 1a Correlation to FACS

  FACS <- read.table("Data/joined_facs_ct", sep="\t", header=T)
    require(plyr)
    func <- function(xx)
      {
        return(data.frame(COR = cor(xx$cell_contribution, xx$percent_of_CD54, use = "complete.obs")))
      }
    FACS_Corr <- ddply(FACS, .(cell_cytoreason, cell_facs), func)
    
    FACS_Corr$cell_cytoreason <-  factor(FACS_Corr$cell_cytoreason, levels =fit$labels[fit$order])
    FACS_Corr$cell_facs <-  factor(FACS_Corr$cell_facs, levels =fit1$labels[fit1$order])

      #Plot heatmap correlation
      A <- ggplot(FACS_Corr, aes(x=cell_cytoreason, y=cell_facs, fill=COR))+geom_tile()+
            scale_colour_gradient2(
              low = "blue",
              mid = "white",
              high = "red",
              midpoint = 0,
              space = "Lab",
              na.value = "grey50",
              guide = "colourbar",
              aesthetics = "fill"
              ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
            geom_rect(mapping=aes(xmin=3.5, xmax=7.5, ymin=0.5, ymax=4.5, fill=NULL), color="black", alpha=0.0)+
            geom_rect(mapping=aes(xmin=7.5, xmax=10.5, ymin=5.5, ymax=7.5, fill=NULL), color="black", alpha=0.0)+
            geom_rect(mapping=aes(xmin=10.5, xmax=15.5, ymin=7.5, ymax=10.5, fill=NULL), color="black", alpha=0.0)+
            geom_rect(mapping=aes(xmin=15.5, xmax=18.5, ymin=12.5, ymax=16.5, fill=NULL), color="black", alpha=0.0)+
            geom_rect(mapping=aes(xmin=18.5, xmax=23.5, ymin=10.5, ymax=12.5, fill=NULL), color="black", alpha=0.0)+
            annotate("text", x=5.5, y=5, label="Monocytes", size=5)+ 
            annotate("text", x=8, y=8, label="NK", size=5)+ 
            annotate("text", x=11.5, y=11, label="CD8", size=5)+ 
            annotate("text", x=17, y=17, label="B Cell", size=5) +
            annotate("text", x=19.5, y=13, label="CD4", size=5)+
            coord_fixed(clip = 'off')

 #Figure 1b
    CellCorPlot <- function(CellDF1 = "Monocyte"){
      CellDF <- get (CellDF1)
      CorLab <- signif(cor(CellDF$cell_contribution, CellDF$percent_of_CD54, use = "complete.obs"),2)
      p <- ggplot(CellDF, aes(y=cell_contribution, x=percent_of_CD54))+geom_point()+ ggtitle(paste0(CellDF1, " ", CorLab))+
        xlab("FACS")+ ylab("Cell Contribution")+theme_cowplot()
      return(p)
          }
    
      Monocyte <- FACS[(FACS$cell_cytoreason=="Monocyte" & FACS$cell_facs=="Monocyte"),]
      matureNK <- FACS[(FACS$cell_cytoreason=="Mature NK" & FACS$cell_facs=="Mature NK"),]
      CD8 <- FACS[(FACS$cell_cytoreason=="CD8+ ab T cell" & FACS$cell_facs=="CD8+ ab T cell"),]
      CD4 <- FACS[(FACS$cell_cytoreason=="Naive CD4+ T cell" & FACS$cell_facs=="Naive CD4+ T cell"),]
      matureBcell <- FACS[(FACS$cell_cytoreason=="Mature B cell" & FACS$cell_facs=="Mature B cell"),]
    
        B4 <- CellCorPlot(CellDF1 = "Monocyte")
        B5 <- CellCorPlot(CellDF1 = "CD8")
        B3 <- CellCorPlot(CellDF1 = "matureNK")
        B2 <- CellCorPlot(CellDF1 = "CD4") 
        B1 <- CellCorPlot(CellDF1 = "matureBcell")   
    
          ComB <- plot_grid(B1, B2, B3, B4, B5, ncol = 2)  
          plot_grid(A, ComB, rel_widths = c(2,1), labels = c('A', 'B'))
          ggsave("Figure1.pdf",  units = "in",  width = 4.5, height = 4.5/3, dpi=300,scale =4)
      
      

## Figure 2 Cell Changes Associated with Disease

  Fig2Df <- rbind(ctDiffMTX[ctDiffMTX$model_term=="adjusted_RA_vs_HC",], 
                  ctDiffTNF[ctDiffTNF$model_term=="adjusted_RA_vs_HC",])
  Fig2Df$Level2 <-c(rep(c(rep("B Cells", 5), rep("NK Cells", 3), rep("T Cells", 14), rep("DC Cells", 3), rep("Monocytes", 4)),2))
   Fig2Df$cell_type <- with(Fig2Df, 
                          gsub(" cell|thymus-derived ", "", 
                               gsub("alpha-beta ", "ab", 
                                    gsub("negative", "-", 
                                        gsub("natural killer", "NK", 
                                            gsub("class switched", "CS",
                                              gsub("-positive", "+", cell_type)))))))


  Levels <- c( "plasma cell", "memory B cell" , "class switched memory B cell", "mature B cell" , "naive B cell",
                    "CD16-positive, CD56-dim natural killer cell", "natural killer cell", "mature natural killer cell",
                    "T-helper 17 cell" ,  "T-helper 1 cell", "central memory CD8-positive, alpha-beta T cell", "CD4-positive, alpha-beta T cell",                            
                    "CD8-positive, alpha-beta T cell",           
                    "effector memory CD4-positive, alpha-beta T cell",            "effector memory CD8-positive, alpha-beta T cell",           
                    "effector memory RA CD8-positive, alpha-beta T cell (TEMRA)" ,                                           
                    "mature NK T cell" ,                     
                    "naive thymus-derived CD4-positive, alpha-beta T cell",       "naive thymus-derived CD8-positive, alpha-beta T cell",      
                    "regulatory T cell" ,                    
                    "T-helper 2 cell",
                    "central memory CD4-positive, alpha-beta T cell", 
                    "myeloid dendritic cell", "dendritic cell" ,  "plasmacytoid dendritic cell" ,
                    "CD14-positive, CD16-negative classical monocyte", "monocyte", "CD14-low, CD16-positive monocyte", "CD14-positive, CD16-positive monocyte"                                     
)

  Fig2Df$cell_type <- factor(Fig2Df$cell_type, levels=unique(Fig2Df$cell_type))

  # New facet label names for dose variable
  dose.labs <- c("Early RA", "Established RA")
  names(dose.labs) <- c("MTX", "TNF")

  ggplot(Fig2Df, aes(x=cell_type, y=-log10(FDR), fill=estimate<0, group="Treatment", label=""))+geom_bar(stat="identity",position = position_dodge2(width=1))+
    facet_wrap(~Treatment,nrow=1, labeller = labeller(Treatment = dose.labs ))+
    coord_flip(clip = "off", ylim = c(0, 2.75))+
    scale_x_discrete(limits = rev(levels(as.factor(Fig2Df$cell_type))))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept = 1., alpha=0.4)

  ggsave("Figure2.png", width = 1.5, height = 1.5,  units = "in", dpi=300,scale =4)
  ggsave("Figure2.pdf", width = 3.0, height = 1.5,  units = "in", dpi=300,scale =4)
    #####
    #Figure 2b Meta Analysis
    # Data downloaded from cytoreason portal on July 9 2020
    #
    Meta2b <- read.table("Data/BloodRAMetaAnalysis.txt", header=T, sep="\t")
      colnames(Meta2b) <- c("cell_type", "FDR", "estimate")
      #remove whole blood specific cells
    #Meta2b$cell_type <- factor(Meta2b$cell_type, levels=Levels)
      WBCells <- c("neutrophil", "granulocyte", "eosinophil")
      Meta3b <- Meta2b[!Meta2b$cell_type%in% WBCells,]
      Meta3b$Treatment ="Meta_analysis"
      Meta3b$cell_type <- with(Meta3b, 
                               gsub(" cell|thymus-derived ", "", 
                                    gsub("alpha-beta ", "ab", 
                                         gsub("negative", "-", 
                                              gsub("natural killer", "NK", 
                                                   gsub("class switched", "CS",
                                                        gsub("-positive", "+", cell_type)))))))
      Meta3b$cell_type <- factor(Meta3b$cell_type, levels =(c("CS memory B" ,"mature B",  "memory B",  "naive B", "plasma",
                                                                "CD16+, CD56-dim NK", "mature NK", "NK",
                                                                "CD4+, abT", "CD8+, abT", "central memory CD4+, abT", "central memory CD8+, abT", "effector memory CD4+, abT", "effector memory CD8+, abT" ,
                                                                "effector memory RA CD8+, abT (TEMRA)",   "mature NK T", "naive CD4+, abT", "naive CD8+, abT", "regulatory T" , 
                                                                "T-helper 1", "T-helper 2", "T-helper 17", "dendritic", "myeloid dendritic"  ,"plasmacytoid dendritic",
                                                                "CD14-low, CD16+ monocyte", "CD14+, CD16-- classical monocyte",   "CD14+, CD16+ monocyte" ,    "monocyte" )))
      
                                                                                                 
                                                                                                                        
                                                                
      
    
      
      ggplot(Meta3b, aes(x=cell_type, y=FDR, fill=estimate<0, group="Treatment", label=""))+
        geom_bar(stat="identity",position = position_dodge2(width=1))+
        facet_wrap(~Treatment,nrow=1)+
        coord_flip(clip = "off", ylim = c(0, 10))+
        scale_x_discrete(limits = rev(levels(as.factor(Meta3b$cell_type))))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept = 1., alpha=0.4)+
        theme(axis.text.y=element_blank())+ylab("")
      
        ggsave("Figure2b.png", width = .66, height = 1.5,  units = "in", dpi=300,scale =4)  
  

####Combined Figure 2
figure2Combined <- function(){
  Meta3b$FDR <- (10^-Meta3b$FDR)
  CombinedFigureData1 <- rbind(Fig2Df[,c(3,8,4,9)], Meta3b)
  # New facet label names for dose variable
  dose.labs <- c("Early RA", "Established RA", "MetaAnalysis")
  names(dose.labs) <- c("MTX", "TNF", "Meta_analysis")
  CombinedFigureData1$Treatment <- factor(CombinedFigureData1$Treatment, levels=c("MTX", "TNF", "Meta_analysis"))

    CombinedFigureData1$cell_type <- with(CombinedFigureData1, 
                                          gsub(" cell|thymus-derived ", "", 
                                               gsub("alpha-beta ", "ab", 
                                                    gsub("negative", "-", 
                                                         gsub("natural killer", "NK", 
                                                              gsub("class switched", "CS",
                                                                   gsub("-positive", "+", cell_type)))))))
    
    
    Level2 <- data.frame(cell_type = Fig2Df$cell_type, Level2 = Fig2Df$Level2)[1:29,]
    Level2$cell_type <- with(Level2, 
                                          gsub(" cell|thymus-derived ", "", 
                                               gsub("alpha-beta ", "ab", 
                                                    gsub("negative", "-", 
                                                         gsub("natural killer", "NK", 
                                                              gsub("class switched", "CS",
                                                                   gsub("-positive", "+", cell_type)))))))
    
    CombinedFigureData1 <- merge(CombinedFigureData1, Level2, by="cell_type")
    
 
  CombinedFigureData1$cell_type <- factor(CombinedFigureData1$cell_type, levels =(c("CS memory B" ,"mature B",  "memory B",  "naive B", "plasma",
                                                                                  "CD16+, CD56-dim NK", "mature NK", "NK",
                                                                                  "CD4+, abT", "CD8+, abT", "central memory CD4+, abT", "central memory CD8+, abT", "effector memory CD4+, abT", "effector memory CD8+, abT" ,
                                                                                  "effector memory RA CD8+, abT (TEMRA)",   "mature NK T", "naive CD4+, abT", "naive CD8+, abT", "regulatory T" , 
                                                                                  "T-helper 1", "T-helper 2", "T-helper 17", "dendritic", "myeloid dendritic"  ,"plasmacytoid dendritic",
                                                                                  "CD14-low, CD16+ monocyte", "CD14+, CD16-- classical monocyte",   "CD14+, CD16+ monocyte" ,    "monocyte" )))


  P <- ggplot(CombinedFigureData1, aes(x=cell_type, y=-log10(FDR), fill=estimate<0, group="Treatment", label=""))+geom_bar(stat="identity",position = position_dodge2(width=1))+
    facet_wrap(~Treatment,nrow=1, labeller = labeller(Treatment = dose.labs ), scales="free_x")+
    #coord_flip(clip = "off", ylim = c(0, 2.75))+
    coord_flip(clip = "off")+
    #scale_x_discrete(limits = rev(levels(as.factor(CombinedFigureData$cell_type))))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept = 1., alpha=0.4)
    # Modify legend titles
    # Edit legend title and labels
  P <- P + scale_fill_discrete(name = "Estimate", labels = c("Increase", "Decrease"))
  print(P)
  #ggsave("Figure2Comb.png", width = 1.5, height = 1.5,  units = "in", dpi=300,scale =4)  
}

figure2Combined()
ggsave("Figure2Comb1.png", width = 1.5, height = 1.5,  units = "in", dpi=300,scale =4)  
ggsave("Figure2Comb1.pdf", width = 1.5, height = 1.5,  units = "in", dpi=300,scale =4)  
###### Clinical Correlation  
######
  ClinCorrfunction <- function(Cell="monocyte", Clin = "Prednisolone_signature"){
    require(cowplot)
    INDEX <- which (rownames(Cellgset)==Cell)
    CorrData <- data.frame(pData(Cellgset)[Clin], exprs(Cellgset)[INDEX,])
      colnames(CorrData) <-c("pData", "exprs")
    CorLab <- signif(cor(CorrData$pData, CorrData$exprs, use = "complete.obs"),2)
    
      p <- ggplot(CorrData, aes(x=pData, y=exprs) )+geom_point() 
        p <- p + theme_half_open()
        p <- p + xlab(Clin) +ylab(Cell) +ggtitle(CorLab)
   return(p)
  }
  a <- ClinCorrfunction(Cell="monocyte", Clin = "Prednisolone_signature")
  b <- ClinCorrfunction(Cell="naive thymus-derived CD8-positive, alpha-beta T cell", Clin = "age")
  plot_grid(a,b) 
  ggsave("Figures/Figure3ClinCorr.pdf", width = 3.0, height = 1.5,  units = "in", dpi=300,scale =4)  
####### End Clinical Correlation  
  
####Figure 4
  
  Figure3Plot <- function(DF = ctDiffMTX, Responders=TRUE){
    DF$model_term <- gsub("adjusted_", "", DF$model_term)
    DF$model_term <- factor(DF$model_term, levels=c("RA_vs_HC", "Post_vs_Pre",
                                                    "Post_vs_Pre_R",   "Post_vs_Pre_NR",
                                                    "baseline_R-NR", "Post_vs_Pre_R-NR"))
   
      DF$cell_type <- gsub(" cell", "", DF$cell_type)
      DF$cell_type <- gsub("-positive", "+",DF$cell_type)
      cellorder <- gsub(" cell", "",cellorder)
      cellorder <- gsub("-positive", "+",cellorder)
      DF$cell_type <- factor(DF$cell_type, levels=cellorder)
    DF <- DF[DF$model_term%in% c("RA_vs_HC", "Post_vs_Pre",
                                "Post_vs_Pre_R",   "Post_vs_Pre_NR"),] 
    DF$Sign <- DF$FDR<0.1
    if (Responders!="TRUE") {
      DF <- DF[DF$model_term%in% c("RA_vs_HC", "Post_vs_Pre"),]
    }
  p3 <- ggplot(DF, aes(x=model_term, y=cell_type, fill=estimate))+geom_tile()+
    scale_colour_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      space = "Lab",
      na.value = "grey50",
      guide = "colourbar",
      aesthetics = "fill"
    )
  p3 <- p3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    geom_point(aes(size=Sign)) +
    scale_size_manual(values=c(-50,3),guide="none")
  #p3 <- p3 +ggtitle("MTX")
  p3  

  

return(p3)
  }

  A <- Figure3Plot(DF = ctDiffMTX)
  A <- A+ggtitle("MTX")
  #ggsave("Figure3c.png", width = 1.5, height = 1.5,  units = "in", dpi=300,scale =4)  
  B <- Figure3Plot(DF = ctDiffTNF)
  B <- B +ggtitle("TNF")
  plot_grid(A,B)
  
  ggsave("Figure4.pdf", width = 4.5, height = 2.5,  units = "in", dpi=300,scale =4) 
  ##### End Figure 3


  
  ####Figure 5
  
  Fig5Fun <- function(DF=gtDiffMTX){
    require(cowplot)
    require(ggrepel)
    library(reshape2)
    Title <- deparse(substitute(DF))
      Title <- gsub("gtDiff", "", Title)
    DF <- DF[DF$term=="Post_vs_Pre",]
      DF <- DF[DF$adjustment_terms %in% unique(DF$adjustment_terms)[c(1,5)],]
      DF.cast <- dcast(DF, SYMBOL ~ adjustment_terms,mean, value.var = "FDR")
        colnames(DF.cast)<- c("SYMBOL", "Adjusted_FDR", "Unadjusted_FDR")
      DF.cast$Color <- ifelse(DF.cast$Adjusted_FDR <0.1 & DF.cast$Unadjusted_FDR <0.1, 1, 
                          ifelse(DF.cast$Adjusted_FDR <0.1 & DF.cast$Unadjusted_FDR >0.1, 2, 
                                  ifelse (DF.cast$Adjusted_FDR >0.1 & DF.cast$Unadjusted_FDR <0.1, 3,4)))
      DF.cast$Alpha <-gsub("1|2|3", 1, DF.cast$Color) 
        DF.cast$Alpha <-gsub("4", 0.2, DF.cast$Alpha) 
        DF.cast$Label <- ifelse(DF.cast$Adjusted_FDR<0.1|DF.cast$Unadjusted_FDR<0.01, as.character(DF.cast$SYMBOL), "")
        
      
    p <- ggplot(DF.cast, aes(x=-log10(Unadjusted_FDR), y=-log10(Adjusted_FDR), alpha = Alpha, color=as.character(Color), label=Label))+geom_point()+ theme_cowplot(12)
    p <- p+ geom_hline(yintercept =1, alpha=0.4, linetype="dashed")+ geom_vline(xintercept =1, alpha=0.4, linetype="dashed")
    p <- p + scale_color_manual(breaks = c("1", "3","2" ,"4"),values=c("red", "blue", "darkgoldenrod", "black"), 
                                labels = c("Before and After Adjustment", "Before Adjustment", "After Adjustment", "Not Sign"))
    p <- p + geom_text_repel(data = subset(DF.cast, Label != ""), show.legend = FALSE)+ ggtitle(Title)
    # Modify legend titles
    p <- p + labs(color = "Significance")
    p <- p + guides(color = FALSE)
    p <- p + guides(alpha = FALSE, label = FALSE)
    p
    return(p)
  }
  A=Fig5Fun(DF=gtDiffMTX)
  B=Fig5Fun(DF=gtDiffTNF)
  plot_grid(A,B)
  
  ggsave("Figure5.pdf", width = 3.0, height = 1.5,  units = "in", dpi=300,scale =4)  
  
  #####End Figure 5
  
 #Table1

TB1 <- read.table("Data/Table1", header=T, sep="\t", as.is = T)

