##
library(Matrix)
library(tidyverse)
library(Seurat)

library(ComplexHeatmap, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(circlize, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(cowplot)
library(ggrepel)
library(ggtext, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(glue)

##

rm(list=ls())

outdir <- "./2_markers.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)



####
###
sc <- read_rds("../3_cluster.outs/1.2_seurat.annot.rds")
meta <- sc@meta.data

### color version 1, "#b3de69", "#fdb462",  "#bc80bd", 
col_MCls <- c("ODC"="#fb8072", "Astrocyte"="#92cd2d", "Microglia"="#c38dc4", "OPC"="#35978f",
   "Non_DA"="#fc9016",  "Pericyte"="#bebada", "DA"="#80b1d3", "Endothelial"="#fccde5",
    "T-cell"="#8dd3c7",  "Ependymal"="#828282")


###
MCls_name <- c("0"="ODC", "1"="Astrocyte", "2"="Microglia", "3"="OPC",
     "4"="Non_DA", "6"="Pericyte", "7"="DA", "8"="Endothelial", "9"="T-cell", "12"="Ependymal")

MCls_val <- c("ODC"=1, "OPC"=2, "Astrocyte"=3, "Microglia"=4, "DA"=5, "Non_DA"=6, 
    "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)



####################################################
### Single cell level for opioid receptors genes ###
####################################################

 
geneList <- c("OPRM1", "OPRK1", "OPRD1", "OPRL1")
plotData <- FetchData(sc,vars=geneList)

plotData<- plotData%>%mutate(MCls=meta$celltype, MCls_value=MCls_val[as.character(MCls)],
                  MCl2=fct_reorder(MCls, as.numeric(MCls_value)))

###
### plot data
figs_ls <- lapply(1:4, function(i){
    ##
    gene <- geneList[i]
    plotDF2 <- plotData[,c(i, 5:7)]
    names(plotDF2)[1] <- "y"
    ##
    p <- ggplot(plotDF2, aes(x=MCl2, y=y, color=MCl2))+
       geom_violin(width=0.5)+
       geom_jitter(width=0.2, size=0.1)+
       ylab("Normalized gene expression")+
       scale_color_manual(values=col_MCls)+ 
       ggtitle(gene)+
       theme_bw()+
       theme(plot.title=element_text(hjust=0.5, size=12),
             axis.text.x=element_text(size=9, angle=45, hjust=1),
             axis.text.y=element_text(size=9),
             axis.title.x=element_blank(),
             axis.title.y=element_text(size=10),
             legend.position="none")
     p
})

### output figures
figfn <- paste(outdir, "Figure1_opioid_violin.png", sep="")
png(figfn, width=680, height=700, res=120)
plot_grid(plotlist=figs_ls, ncol=2)
dev.off()




###########################################
### bulk data for opioid receptor genes ###
###########################################

fn <- "../4.0_Diff.cluster.outs/YtX_sel.comb.rds"
x <- read_rds(fn)
depth <- colSums(x)
x2 <- sweep(x, 2, depth, "/")*1e+6

##
cvt <- str_split(colnames(x2), "_", simplify=T)
plotDF <- data.frame(rn=colnames(x2), MCls=MCls_name[cvt[,1]], sampleID=cvt[,2]) 
### 
plotDF <- plotDF%>%mutate(MCls_value=MCls_val[as.character(MCls)],
                          MCl2=fct_reorder(MCls, as.numeric(MCls_value)))

###
figs_ls <- lapply(geneList, function(gene){
    ##
    cat(gene, "\n")
    plotDF2 <- plotDF%>%mutate(y=log2(x2[gene,]+1))
    ###
    p <- ggplot(plotDF2, aes(x=MCl2, y=y, color=MCl2))+
       geom_violin(width=0.5)+
       geom_jitter(width=0.2, size=0.1)+
       ylab("Normalized gene expression")+
       scale_color_manual(values=col_MCls)+ 
       ggtitle(gene)+
       theme_bw()+
       theme(plot.title=element_text(hjust=0.5, size=12),
             axis.text.x=element_text(size=9, angle=45, hjust=1),
             axis.text.y=element_text(size=9),
             axis.title.x=element_blank(),
             axis.title.y=element_text(size=10),
             legend.position="none")
     p
})

### output figures
figfn <- paste(outdir, "Figure1.2_opioid_indi_violin.png", sep="")
png(figfn, width=680, height=700, res=120)
plot_grid(plotlist=figs_ls, ncol=2)
dev.off()


## DF <- map_dfr(geneList, function(ii){
##    ##
##    plotDF2 <- plotDF%>%mutate(y=x2[ii,], y0=log2(y), y2=log2(y+1), gene=ii)
##    plotDF2
## })


