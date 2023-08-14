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

rm(list=ls())


##
outdir <- "./2_seurat.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###
### fun, read cellranger outputs into seurat object
readCellRanger <- function(basefolder, exp){
###    
   countfn <- paste(basefolder, exp, "/outs/filtered_feature_bc_matrix", sep="")
   counts <- Read10X(countfn)
   sc <- CreateSeuratObject(counts, project="sc-brain")
   newcell <- paste(exp, "_", gsub("-1", "", Cells(sc)), sep="")
   sc <- RenameCells(sc, new.names=newcell)
   sc
}



##############################################
### create seurat object for 6 experiments ###
##############################################


basefolder <- "/wsu/home/groups/bannonlab/sc_brains/counts_brain_2022_merge/"
fn <- paste(basefolder, "libList.txt", sep="")
libList <- read.table(fn)$V1

sc_ls <- mclapply(libList, function(ii){
    cat(ii, "\n")
    sc <- readCellRanger(basefolder, ii)
    sc
},mc.cores=16)

combined <- merge(sc_ls[[1]], sc_ls[-1], project="sc-brain")
x <- combined@meta.data
x$NEW_BARCODE <- rownames(x)
combined <- AddMetaData(combined, x)
###
opfn <- paste(outdir, "1.0_seurat.all.rds", sep="")
write_rds(combined, opfn)



### add demuxlet results
fn <- paste(outdir, "1.0_seurat.all.rds", sep="")
sc <- read_rds(fn)
x <- sc@meta.data
###
meta <- read_rds("./1_demux_output/1_demux_New.ALL.rds")
meta2 <- meta%>%dplyr::select("DROPLET.TYPE", "SNG.BEST.GUESS", "NEW_BARCODE", "EXP")
x <- x%>%left_join(meta2, by="NEW_BARCODE")
rownames(x) <- x$NEW_BARCODE

combined <- AddMetaData(combined, x)

opfn <- paste(outdir, "1.0_seurat.all.rds", sep="")
write_rds(combined, opfn)


### SNG
fn <- paste(outdir, "1.0_seurat.all.rds", sep="")
sc <- read_rds(fn)
sc2 <- subset(sc, subset=DROPLET.TYPE=="SNG")
###
sc2[["percent.mt"]] <- PercentageFeatureSet(sc2, pattern="^MT-")

opfn <- paste(outdir, "1.1_seurat.SNG.rds", sep="")
write_rds(sc2, opfn)

    
##################################################
### summary of quality data before clean data  ###
##################################################

fn <- paste(outdir, "1.1_seurat.SNG.rds", sep="")
sc <- read_rds(fn)
x <- sc@meta.data

## ### droplet type
## dd <- x%>%
##    group_by(EXP, DROPLET.TYPE)%>%
##    summarise(ncell=n(),.groups="drop")

## dd <- dd%>%drop_na(EXP) ## I will check NA values

## d2 <- dd%>%group_by(EXP)%>%mutate(perc=ncell/sum(ncell))%>%ungroup()
## d2%>%filter(DROPLET.TYPE=="DBL")


## p <- ggplot(dd)+
##    geom_bar(stat="identity", position=position_fill(reverse=F),
##         aes(x=EXP, y=ncell, fill=DROPLET.TYPE))+
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          axis.title=element_blank(),
##          axis.text.x=element_text(angle=90, hjust=1, size=8))

## png("./2_seurat.outs/Figure0.1_droplet.png", width=420, height=420, res=120)
## p
## dev.off()


### number of cells
dd <- x%>%group_by(sampleID)%>%
   summarise(ncell=n(), nreads=mean(nCount_RNA), ngenes=mean(nFeature_RNA), .groups="drop")

EXP_df <- data.frame(EXP=sort(unique(x$EXP)))%>%
   mutate(values=as.numeric(gsub("B|M|R", "", EXP)))
EXP_df <- EXP_df%>%mutate(value2=ifelse(grepl("^M", EXP), values+12, values))

dd <- dd%>%mutate(Batch=str_sub(EXP, 1, 1))%>%
    left_join(EXP_df, by="EXP")%>%
    mutate(EXP=fct_reorder(EXP, value2))%>%as.data.frame()

p <- ggplot(dd)+
    geom_bar(stat="identity", aes(x=EXP, y=ncell, fill=factor(Batch)))+
    ggtitle("Number of cells")+
    theme_bw()+
    theme(legend.position="none",
       axis.title=element_blank(),
       axis.text.x=element_text(angle=90, hjust=1, size=8),
       plot.title=element_text(hjust=0.5))

png("./2_seurat.outs/Figure0.2_cells.png", width=550, height=420, res=120)
p
dev.off()


### reads per cell
p <- ggplot(dd)+
    geom_bar(stat="identity", aes(x=EXP, y=nreads, fill=factor(Batch)))+
    ggtitle("Reads per cell")+
    theme_bw()+
    theme(legend.position="none",
       axis.title=element_blank(),
       axis.text.x=element_text(angle=90, hjust=1, size=8),
       plot.title=element_text(hjust=0.5))

png("./2_seurat.outs/Figure0.3_reads.png", width=550, height=420, res=120)
p
dev.off()

###
### genes per cell
p <- ggplot(dd)+
    geom_bar(stat="identity", aes(x=EXP, y=ngenes, fill=factor(Batch)))+
    ggtitle("Genes per cell")+
    theme_bw()+
    theme(legend.position="none",
       axis.title=element_blank(),
       axis.text.x=element_text(angle=90, hjust=1, size=8),
       plot.title=element_text(hjust=0.5))

png("./2_seurat.outs/Figure0.4_genes.png", width=550, height=420, res=120)
p
dev.off()


###
### violin plots show percent.mt
dd2 <- x%>%left_join(EXP_df, by="EXP")%>%
   mutate(EXP=gsub("R", "", EXP), EXP=fct_reorder(EXP, value2), Batch=str_sub(EXP, 1, 1))

p2 <- ggplot(dd2, aes(x=EXP, y=percent.mt, fill=factor(Batch)))+
    geom_violin(width=1)+
    ## geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
    ylab("percentage of mitochondrial counts")+
    ## facet_wrap(~factor(Batch), nrow=2, labeller=as_labeller(c("B"="Bannon", "M"="Mash")))+
    theme_bw()+
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_text(size=8))
###
figfn <- paste(outdir, "Figure0.5_mt.violin.png", sep="")
png(figfn, width=850, height=400, res=120)
print(p2)
dev.off()


###
###
## dd <- x%>%group_by(SNG.BEST.GUESS, EXP)%>%summarise(ncell=n(), .groups="drop")
## dd2 <- dd%>%group_by(SNG.BEST.GUESS)%>%
##    mutate(Perc=ncell/sum(ncell))

## p <- ggplot(dd2, aes(x=SNG.BEST.GUESS, y=EXP, fill=Perc))+
##    geom_tile()+
##    scale_fill_gradient("Fraction of cells",
##       low="#ffffc8", high="#7d0025", na.value=NA)+
##    theme_bw()+
##    theme(axis.text.x=element_text(hjust=1, vjust=0.5, angle=90, size=6),
##          axis.text.y=element_text(size=8),
##          axis.title=element_blank())

## ###
## figfn <- "./2_seurat.outs/Figure0.5_heatmap.png"
## png(figfn, width=600, height=280, res=100)
## print(p)
## dev.off()

## x2 <- x%>%select(sampleID=SNG.BEST.GUESS, EXP)
## dd <- x2%>%group_by(sampleID, EXP)%>%
##    summarise(ncell=n(), .groups="drop")


## ###
## dd2 <- dd%>%group_by(sampleID)%>%top_n(1, wt=ncell)%>%
##     mutate(comb=paste(EXP, sampleID, sep="_"))
## dd2 <- dd2%>%as.data.frame()%>%dplyr::select(sampleID, comb)

## dd <- dd%>%left_join(dd2, by="sampleID")
## opfn <- "./2_seurat.outs/sample_EXP.comb.rds"
## write_rds(dd, opfn)

 
## ## dd <- dd%>%drop_na(EXP)

## ### data for heatmap
## dd3 <- dd%>%group_by(sampleID)%>%
##    mutate(Perc=ncell/sum(ncell))%>%as.data.frame()

## ### x labels
## ## sampleID <- dd2$sampleID
## ## names(sampleID) <- dd2$comb
## p <- ggplot(dd3, aes(x=comb, y=EXP, fill=Perc))+
##    geom_tile()+
##    scale_fill_gradient("Fraction of cells",
##       low="#ffffc8", high="#7d0025", na.value=NA)+
##    theme_bw()+
##    theme(axis.text.x=element_blank(), #element_text(hjust=1, vjust=0.5, angle=90, size=6),
##          axis.text.y=element_text(size=8),
##          axis.title=element_blank())

## ###
## figfn <- paste(outdir, "Figure2.5_heatmap.png", sep="")
## png(figfn, width=600, height=300, res=100)
## print(p)
## dev.off()



##################
### clean data ###
##################

## remove the droplets, (1). individual not match library, (2). mt>=20.

sc <- read_rds("./2_seurat.outs/1.1_seurat.SNG.rds")

x <- sc@meta.data
x <- x%>%mutate(sampleID=toupper(gsub(".*_", "", SNG.BEST.GUESS)),
                comb=paste(sampleID, EXP, sep="_"),
                Batch=ifelse(grepl("AKB", sampleID), "Bannon", "Mash"))
sc <- AddMetaData(sc, x)

## correct sample-library
brainCV <- openxlsx::read.xlsx("cell.counts_2022-07-06_final.xlsx")
brainCV <- brainCV[,1:17]
##
brainCV <- brainCV%>%filter(USE==1)%>%
    mutate(Sample_ID2=toupper(Sample_ID2), comb=paste(Sample_ID2, Library, sep="_"))

## subset
cellSel <- x%>%filter(comb%in%unique(brainCV$comb), percent.mt<20)%>%pull(NEW_BARCODE)
sc2 <- subset(sc, cells=cellSel) ## clean data, 36601 genes*212,713 cells

##
opfn <- "./2_seurat.outs/2_seurat.clean.rds"
write_rds(sc2, opfn)



###
### summary clean data
## sc <- read_rds("./2_seurat.outs/2_seurat.clean.rds")
## x <- sc@meta.data

### number of cells
## dd <- x%>%group_by(EXP)%>%
##    summarise(ncell=n(), nreads=mean(nCount_RNA), ngenes=mean(nFeature_RNA), .groups="drop")

## EXP_df <- data.frame(EXP=sort(unique(x$EXP)))%>%
##    mutate(values=as.numeric(gsub("B|M|R", "", EXP)))
## EXP_df <- EXP_df%>%mutate(value2=ifelse(grepl("^M", EXP), values+12, values))

## dd <- dd%>%mutate(Batch=str_sub(EXP, 1, 1))%>%
##     left_join(EXP_df, by="EXP")%>%
##     mutate(EXP=fct_reorder(EXP, value2))%>%as.data.frame()

## ##
## ### number of individual per library
## dx <- x%>%distinct(sampleID, .keep_all=T)
## dx <- dx%>%group_by(EXP)%>%summarise(nindi=n(),.groups="drop")%>%
##    mutate(Batch=str_sub(EXP, 1, 1))%>%
##    as.data.frame()
## ##
## dx <- dx%>%left_join(EXP_df, by="EXP")%>%
##     mutate(EXP=fct_reorder(EXP, value2))

## p1 <- ggplot(dx)+
##     geom_bar(stat="identity",  aes(x=EXP, y=nindi, fill=factor(Batch)))+
##     ylab("#individual per library")+
##     theme_bw()+
##     theme(legend.position="none",
##           axis.title.x=element_blank(),
##           axis.text.x=element_text(angle=90, hjust=1, size=8))
## ##
## png("./2_seurat.outs/Figure1.1_indi.png", width=550, height=420, res=120)
## p1
## dev.off()
          



## p <- ggplot(dd)+
##     geom_bar(stat="identity", aes(x=EXP, y=ncell, fill=factor(Batch)))+
##     ggtitle("Number of cells")+
##     theme_bw()+
##     theme(legend.position="none",
##        axis.title=element_blank(),
##        axis.text.x=element_text(angle=90, hjust=1, size=8),
##        plot.title=element_text(hjust=0.5))

## png("./2_seurat.outs/Figure1.2_cells.png", width=550, height=420, res=120)
## p
## dev.off()


## ### reads per cell
## p <- ggplot(dd)+
##     geom_bar(stat="identity", aes(x=EXP, y=nreads, fill=factor(Batch)))+
##     ggtitle("Reads per cell")+
##     theme_bw()+
##     theme(legend.position="none",
##        axis.title=element_blank(),
##        axis.text.x=element_text(angle=90, hjust=1, size=8),
##        plot.title=element_text(hjust=0.5))

## png("./2_seurat.outs/Figure1.3_reads.png", width=550, height=420, res=120)
## p
## dev.off()

## ###
## ### genes per cell
## p <- ggplot(dd)+
##     geom_bar(stat="identity", aes(x=EXP, y=ngenes, fill=factor(Batch)))+
##     ggtitle("Genes per cell")+
##     theme_bw()+
##     theme(legend.position="none",
##        axis.title=element_blank(),
##        axis.text.x=element_text(angle=90, hjust=1, size=8),
##        plot.title=element_text(hjust=0.5))

## png("./2_seurat.outs/Figure1.4_genes.png", width=550, height=420, res=120)
## p
## dev.off()

## counts <- sc@assays$RNA@counts
## nreads <- colSums(counts)
## ngenes <- colSums(counts>0)

## x2 <- x%>%dplyr::select(NEW_BARCODE, sampleID)%>%
##     mutate(nCount_RNA=colSums(counts), nFeature_RNA=colSums(counts>0))
## dd2 <- x2%>%group_by(sampleID)%>%
##    summarise(ncell=n(), nreads=mean(nCount_RNA), ngenes=mean(nFeature_RNA), .groups="drop")

####################################################
### compare the number of cells, genes and depth ###
####################################################


sc <- read_rds("./3_cluster.outs/1.2_seurat.annot.rds")
x <- sc@meta.data

### number of cells
cvt <- openxlsx::read.xlsx("BrainCV_cell_counts_final_correct_2022-12-06.xlsx")
cvt2 <- cvt%>%dplyr::select(sampleID, Category)
###
dd <- x%>%group_by(sampleID)%>%
   summarise(ncell=n(), nreads=mean(nCount_RNA), ngenes=mean(nFeature_RNA), .groups="drop")
##
dd <- dd%>%left_join(cvt2, by="sampleID")


## cell-type-sampleID
dd2 <- x%>%group_by(sampleID, celltype)%>%summarise(ny=n(), .groups="drop")%>%
    pivot_wider(id_cols=sampleID, names_from=celltype, values_from=ny, values_fill=0)
dd2 <- dd2[,c("sampleID", "ODC", "OPC", "Astrocyte", "Microglia", "DA", "Non_DA",
              "Pericyte", "Endothelial", "T-cell", "Ependymal")]

dd <- dd%>%left_join(dd2, by="sampleID")

###
###
fn <- paste(outdir, "TableS_Data_quality.xlsx", sep="")
openxlsx::write.xlsx(dd, file=fn, overwrite=T)


###
### test the difference between control and opioid
dd1 <- dd%>%dplyr::select(y=ncell, x=Category)
y1 <- dd1%>%filter(x=="Control")%>%pull(y)
y2 <- dd1%>%filter(x=="Opioid")%>%pull(y)
pval <- t.test(y1, y2)$p.value

###
## if (pval>=0.05){
##     symbol <- "ns"
## }else if (pval<0.05 & pval>0.01){
##     symbol <- "*"
## }else if (pval<0.01&p>0.001){
##     symbol <- "**"
## }else{
##     symbol <- "***"
## }
    
p1 <- ggplot(dd1, aes(x=factor(x), y=y, fill=x))+
   geom_violin(width=0.8)+
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
   scale_fill_manual(values=c("Control"="#4575b4", "Opioid"="#d73027"))+ 
   geom_signif(comparison=list(c("Control", "Opioid")),
               annotation="NS",
               y_position=c(max(y1)+100, max(y2)),
               tip_length=0.05, vjust=0, textsize=3)+
   scale_y_continuous("#Cells per sample",
                      expand=expansion(mult=c(0.3, 0.3)))+ 
   theme_bw()+
   theme(axis.title=element_text(size=10),
         axis.title.x=element_blank(),
         axis.text=element_text(size=10),
         ##plot.title=element_text(hjust=0.5, size=12),
         legend.position="none")


###
###  #genes
    
dd2 <- dd%>%dplyr::select(y=ngenes, x=Category)
y1 <- dd2%>%filter(x=="Control")%>%pull(y)
y2 <- dd2%>%filter(x=="Opioid")%>%pull(y)
pval <- t.test(y1, y2)$p.value

##
p2 <- ggplot(dd2, aes(x=factor(x), y=y, fill=x))+
   geom_violin(width=0.8)+
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
   scale_fill_manual(values=c("Control"="#4575b4", "Opioid"="#d73027"))+ 
   geom_signif(comparison=list(c("Control", "Opioid")),
               annotation="NS",
               y_position=c(max(y1)+50, max(y2)),
               tip_length=0.05, vjust=0, textsize=3)+
   scale_y_continuous("Average #genes per cell for each sample",
                      expand=expansion(mult=c(0.3,0.3)))+  
   theme_bw()+
   theme(axis.title.y=element_text(size=10),
         axis.title.x=element_blank(),
         axis.text=element_text(size=10),
         ##plot.title=element_text(hjust=0.5, size=12),
         legend.position="none")

###
### depeth
    
dd3 <- dd%>%dplyr::select(y=nreads, x=Category)
y1 <- dd3%>%filter(x=="Control")%>%pull(y)
y2 <- dd3%>%filter(x=="Opioid")%>%pull(y)
pval <- t.test(y1, y2)$p.value

    
p3 <- ggplot(dd3, aes(x=factor(x), y=y, fill=x))+
   geom_violin(width=0.8)+
   geom_boxplot(width=0.2, color="grey", outlier.shape=NA)+
   scale_fill_manual(values=c("Control"="#4575b4", "Opioid"="#d73027"))+ 
   geom_signif(comparison=list(c("Control", "Opioid")),
               annotation="NS",
               y_position=c(max(y1)+200, max(y2)),
               tip_length=0.05, vjust=0, textsize=3)+
   scale_y_continuous("Average #reads per cell for eachsample",
                      expand=expansion(mult=c(0.3,0.3)))+ 
   theme_bw()+
   theme(axis.title.y=element_text(size=10),
         axis.title.x=element_blank(),
         axis.text=element_text(size=10),
         ##plot.title=element_text(hjust=0.5, size=12),
         legend.position="none")

figfn <- paste(outdir, "Figure1_qc.png", sep="")
png(figfn, width=900, height=380, res=120)    
plot_grid(p1, p2, p3, labels="AUTO", label_fontface="plain", ncol=3, label_x=0.1, label_y=1.01)
dev.off()
    
    



####################################################################
### heatmap show distribution of individual in different batches ###
####################################################################

fn <- paste(outdir, "2_seurat.clean.rds", sep="")
sc <- read_rds(fn)
x <- sc@meta.data

dd <- x%>%group_by(sampleID, EXP)%>%summarise(ncell=n(), .groups="drop")%>%as.data.frame()

dd <- dd%>%mutate(yvalue=as.numeric(gsub("^[BM]|R", "", EXP)),
            yvalue=ifelse(grepl("^M", EXP), yvalue+12, yvalue))

##
cvt <- openxlsx::read.xlsx("BrainCV_cell_counts_final_correct_2022-12-06.xlsx")
cvt2 <- cvt%>%dplyr::select(ExpectedPool, sampleID)%>%    
   mutate(xvalue=as.numeric(gsub("^[BM]|R", "", ExpectedPool)),
           xvalue=ifelse(grepl("^M", ExpectedPool), xvalue+12, xvalue),
          sampleID2=paste(ExpectedPool, sampleID, sep="_"))
    

## brainCV <- read.csv("BrainCV_cell.counts.csv")
## brainCV <- brainCV%>%filter(USE==1)%>%
##     mutate(Sample_ID2=toupper(Sample_ID2), comb=paste(Sample_ID2, Library, sep="_"))
## ###
## brainCV2 <- brainCV%>%
##    dplyr::select(EXP2=ExpectedPool, sampleID=Sample_ID2)%>%
##    mutate(xvalue=as.numeric(gsub("^[BM]|R", "", EXP2)),
##            xvalue=ifelse(grepl("^M", EXP2), xvalue+12, xvalue),
##           sampleID2=paste(EXP2, sampleID, sep="_"))

## brainCV2 <- brainCV%>%dplyr::select(sampleID=Sample_ID2, Library, Category)



dd2 <- dd%>%left_join(cvt2, by="sampleID")
dd2 <- dd2%>%mutate(EXP=fct_reorder(EXP, yvalue), sampleID2=fct_reorder(sampleID2, xvalue))

## dd2 <- dd2%>%mutate(datasets=ifelse(grepl("^B", EXP), "Bannon", "Mash"))
## dd2 <- dd <- dd%>%left_join(brainCV2, by="sampleID")
##



### data for heatmap
## dd3 <- dd2%>%group_by(sampleID)%>%
##    mutate(Perc=ncell/sum(ncell))%>%as.data.frame()

## ### x labels
## ## sampleID <- dd2$sampleID
## ## names(sampleID) <- dd2$comb
## p <- ggplot(dd3, aes(x=sampleID2, y=EXP, fill=Perc))+
##    geom_tile()+
##    scale_fill_gradient("Percents nuclei",
##       low="#ffffcc", high="#e31a1c", na.value=NA, limits=c(0, 0.8),
##       guide=guide_legend(keywidth=grid::unit(0.5,"cm"), keyheight=grid::unit(1, "cm")))+
##    ylab("Observed library")+ 
##    theme_bw()+
##    theme(axis.text.x=element_text(hjust=0, vjust=0.5, angle=-90, size=5),
##          axis.text.y=element_text(size=10),
##          axis.title.x=element_blank(),
##          axis.title.y=element_text(size=12))

## ###
## figfn <- "./2_seurat.outs/Figure2.1_heatmap.png"
## png(figfn, width=1000, height=550, res=120)
## print(p)
## dev.off()


###
dd4 <- dd2%>%group_by(EXP)%>%
   mutate(Perc=ncell/sum(ncell))%>%as.data.frame()

### x labels
## sampleID <- dd2$sampleID
## names(sampleID) <- dd2$comb
p2 <- ggplot(dd4, aes(x=sampleID2, y=EXP, fill=Perc))+
   geom_tile()+
   geom_tile()+
   scale_fill_gradient("Percent nuclei",
      low="#ffffcc", high="#e31a1c", na.value=NA, limits=c(0, 0.8),
      guide=guide_legend(keywidth=grid::unit(0.6,"cm"), keyheight=grid::unit(0.8, "cm")))+
   ylab("Observed library")+ 
   theme_bw()+
   theme(axis.text.x=element_text(hjust=0, vjust=0.5, angle=-90, size=5),
         axis.text.y=element_text(size=10),
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12))

###
figfn <- "./2_seurat.outs/plots_pub/Figure2.2_heatmap.png"
png(figfn, width=1000, height=550, res=120)
print(p2)
dev.off()


###
### reads mapping to X and Y chromosome
##library(annotables)


### End


