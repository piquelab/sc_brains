##
library(Matrix)
library(tidyverse)
library(data.table)
library(Seurat)
library(parallel)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
##
library(ggrastr)
library(ggsignif)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(openxlsx)

rm(list=ls())

 
##
outdir <- "./1_data_quality.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


####################################################################################
### summary clean data, correlation analysis between covariates and data quality ### 
####################################################################################

###
### meta data
sc <- read_rds("../3_cluster.outs/1.2_seurat.annot.rds")
x <- sc@meta.data


###
### covariates
cvt <- openxlsx::read.xlsx("../BrainCV_cell_counts_final_correct_2022-12-06.xlsx")
cvt2 <- cvt%>%dplyr::select(-ExpectedPool, -CauseOfDeath, -BrainRegion, -Match, -USE, -n)


###
### 
dd <- x%>%group_by(sampleID)%>%
   summarise(ncell=n(), nreads=mean(nCount_RNA), ngenes=mean(nFeature_RNA), mt=mean(percent.mt), .groups="drop")

##
dd <- dd%>%left_join(cvt2, by="sampleID")


### function
feq <- function(x){
  r <- round(as.numeric(x$estimate),digits=3)
  p <- x$p.value
  if(p<0.001) symb <- "***"
  if(p>=0.001 & p<0.01) symb <- "**"
  if (p>=0.01 & p<0.05) symb <- "*"
  if(p>0.05) symb <- "NS"
  
  eq <- deparse(bquote(italic(R)==.(r)~","~.(symb))) ##deparse#)
  eq 
}

## ### x position
## xFun <- function(dx, a=0.7){
##    ### 
##    R <- max(dx$x)-min(dx$x)
##    xpos <- min(dx$x)+a*R
##    xpos
## }    

## ## y position
## yFun <- function(dx, a=0.9){
##    ##
##    R <- max(dx$y)-max(dx$y)
##    ypos <- min(dx$y)+a*R
##    ypos
## }    



## #### facet labels
## qcs_lab <- as_labeller(c("ncell"="#Nuclei per individual",
##            "nreads"="#Genes/nuclei for individual",
##            "ngenes"="#Reads/nuclei for individual",
##            "mt"="#mitochondrial genes(%)/nuclei"))


plotDF <- dd%>%dplyr::select(ncell, nreads, ngenes, mt, Age, pH, PMI, PC1, PC2)

cov_lab <- c("Age", "pH", "PMI", "PC1", "PC2")
qc_lab <- c("ncell", "nreads", "ngenes", "mt")


###plot data
figs_ls <- list()
nn <- 0
for ( k in 1:5){
   kk <- k+4
   for (i in 1:4){
      ###
      nn <- nn+1 
      df <- plotDF[,c(i,kk)]
      names(df) <- c("y", "x")

      df <- df%>%drop_na(x, y)
       
      Rx <- max(df$x)-min(df$x) 
      x0 <- min(df$x)+0.7*Rx
      Ry <- max(df$y)-min(df$y)
      y0 <- min(df$y)+1.2*Ry

      corr <- cor.test(df$x, df$y)
      eq <- feq(corr)  
       
      ##
      p0 <- ggplot(df, aes(x=x, y=y))+
         geom_point(size=1, color="grey50")+
         annotate("text", x=x0, y=y0, label=eq, color="red", parse=T, size=3)+
         scale_x_continuous(expand=expansion(mult=c(0, 0.2)))+
         scale_y_continuous(expand=expansion(mult=c(0, 0.4)))+                                              
         theme_bw()+
         theme(axis.title=element_blank(),
               axis.text=element_text(size=10))
       
      figs_ls[[nn]] <- p0
       
      cat(cov_lab[k], qc_lab[i], "\n")
   } ### End for loop qc
} ### End for loop covs  



## figfn <- paste(outdir, "Figure1_qc_comb.png", sep="")
## png(figfn, width=1050, height=800, res=120)
## plot_grid(plotlist=figs_ls, nrow=5, ncol=4)
## dev.off()



plotDF <- dd%>%dplyr::select(ncell, nreads, ngenes, mt, Sex)
nn <- 20
for (i in 1:4){
    ###
    nn <- nn+1
    df <- plotDF[,c(i,5)]
    names(df) <- c("y", "x")
    
    ### test the difference between control and opioid
    y1 <- df%>%filter(x=="F")%>%pull(y)
    y2 <- df%>%filter(x=="M")%>%pull(y)
    pval <- t.test(y1, y2)$p.value
    
   ### label significance
   if (pval>=0.05){
      symbol <- "ns"
   }else if (pval<0.05 & pval>0.01){
      symbol <- "*"
   }else if (pval<0.01&p>0.001){
      symbol <- "**"
   }else{
      symbol <- "***"
   }

   R1 <- 0.02*max(y1) 
   R2 <- 0.02*max(y2)
    
   p0 <- ggplot(df, aes(x=factor(x), y=y, fill=x))+
      geom_violin(width=0.8)+
      geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
      scale_fill_manual(values=c("F"="#f1b6da", "M"="#b8e186"))+ 
      geom_signif(comparison=list(c("F", "M")),
               annotation=symbol,
               y_position=c(max(y1)+R1, max(y2)+R2),
               tip_length=0.05, vjust=0, textsize=3)+
      scale_x_discrete(labels=c("F"="Female", "M"="Male"))+
      scale_y_continuous(expand=expansion(mult=c(0.3,0.3)))+ 
      theme_bw()+
      theme(axis.title=element_blank(),
            axis.text=element_text(size=10),
            legend.position="none")

    figs_ls[[nn]] <- p0
}





####
#### compare data quality with library

plotDF <- dd%>%dplyr::select(ncell, nreads, ngenes, mt, Library)

nn <- 24
for (i in 1:4){
    ###
    nn <- nn+1
    df <- plotDF[,c(i,5)]
    names(df) <- c("y", "x")
    df <- df%>%mutate(gr=ifelse(grepl("B",x), "gr1", "gr2"))
    df <- df%>%group_by(x)%>%
        mutate(x_val=median(y))%>%ungroup()%>%
        mutate(x2=fct_reorder(x, x_val)) 
    
    x0 <- 10
    Ry <- max(df$y)-min(df$y)
    y0 <- min(df$y)+1.2*Ry
    
   p0 <- ggplot(df, aes(x=factor(x2), y=y))+
      geom_jitter(width=0.05, height=0.1, color="grey", size=1)+
      annotate("text", x=x0, y=y0, label="ANOVA ***", color="red", size=3)+
      ## scale_color_manual(values=c("gr1"="grey", "gr2"="grey30"),
      ##                 labels=c("gr1"="Detroit", "gr2"="Miami"))+ 
      scale_y_continuous(expand=expansion(mult=c(0.3,0.3)))+ 
      theme_bw()+
      theme(axis.title=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=10),
            axis.ticks.x=element_blank(), legend.position="none")    
    ##
    figs_ls[[nn]] <- p0
}

## figfn <- paste(outdir, "Figure1.2_qc_comb.png", sep="")
## png(figfn, width=1000, height=300, res=120)
## plot_grid(plotlist=figs_ls,  ncol=4)
## dev.off()


figfn <- paste(outdir, "Figure1_qc_comb.png", sep="")
png(figfn, width=1150, height=1300, res=120)
plot_grid(plotlist=figs_ls, nrow=7, ncol=4)
dev.off()


summ <- NULL
for ( i in 1:4){
   ## 
   df <- plotDF[,c(i, 5)]
   names(df) <- c("y", "x") 
   lm0 <- lm(y~factor(x), data=df)
   aa <- as.data.frame(anova(lm0))
   aa$qc <- qc_lab[i]
   summ <- rbind(summ, aa)
}

###
opfn <- paste(outdir, "2_qc_anova_lib.xlsx", sep="")
write.xlsx(summ, opfn)

    


## plotDF <- plotDF%>%drop_na(x,y)%>%
##     mutate(cov_value=as.numeric(cov_val[as.character(covs)]),
##            qc_value=as.numeric(qc_val[as.character(qcs)]),
##            cov2=fct_reorder(covs, cov_value),
##            qc2=fct_reorder(qcs, qc_value))

## annoDF <- plotDF%>%group_by(cov2, qc2)%>%nest()%>%
##     mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y)),
##            eq=map(corr, feq),
##            xpos=map_dbl(data, ~xFun(.x, a=0.7)),
##            ypos=map_dbl(data, ~yFun(.x, a=0.9)))%>%
##     dplyr::select(-data, -corr)

       
## ##
## p0 <- ggplot(plotDF, aes(x=x, y=y))+
##     geom_point(size=1, color="grey50")+
##     geom_text(data=annoDF, aes(x=xpos, y=ypos,label=eq), color="red", size=2.5, parse=T)+
##     facet_grid(cov2~qc2, scales="free",
##                labeller=labeller(qc2=qcs_lab))+
##     theme_bw()+
##     theme(axis.title=element_blank(),
##           axis.text=element_text(size=9),
##           strip.text=element_text(size=12),
##           strip.background=element_blank())


###
###






##
## remove the droplets, (1). individual not match library

sc <- read_rds("./2_seurat.outs/1.1_seurat.SNG.rds")

x <- sc@meta.data
x <- x%>%mutate(sampleID=toupper(gsub(".*_", "", SNG.BEST.GUESS)),
                comb=paste(sampleID, EXP, sep="_"),
                Batch=ifelse(grepl("AKB", sampleID), "Bannon", "Mash"))
sc <- AddMetaData(sc, x)

## correct sample-library
brainCV <- read.xlsx("BrainCV_cell_counts_final_correct_2022-12-06.xlsx")
##
brainCV <- brainCV%>%filter(USE==1)%>%
    mutate(sampleID2=toupper(sampleID), comb=paste(sampleID2, Library, sep="_"))

## subset 
x2 <- x%>%dplyr::filter(comb%in%unique(brainCV$comb))

opfn <- paste(outdir, "1.1_meta_match.rds", sep="")
write_rds(x2, opfn)



### histogram distribution
x2 <- read_rds("./2_seurat.outs/response_to_reviews/1.1_meta_match.rds")

p1 <- ggplot(x2, aes(x=percent.mt))+
   geom_histogram(color="grey", fill=NA)+
   geom_vline(xintercept=c(10, 20), linetype="dashed", bins=100, size=0.25, color="red")+
   xlab("percent of mitochondrial gene")+
   ylab("#Nuclei")+
   theme_bw()+
   theme(axis.text=element_text(size=12),
         axis.title=element_text(size=12))

figfn <- paste(outdir, "Figure2_mt_hist.png", sep="")
png(figfn, width=420, height=320, res=120)
print(p1)
dev.off()




