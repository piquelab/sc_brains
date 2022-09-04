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
vx$NEW_BARCODE <- rownames(x)
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
dd <- x%>%group_by(EXP)%>%
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
sc <- read_rds("./2_seurat.outs/2_seurat.clean.rds")
x <- sc@meta.data

### number of cells
dd <- x%>%group_by(EXP)%>%
   summarise(ncell=n(), nreads=mean(nCount_RNA), ngenes=mean(nFeature_RNA), .groups="drop")

EXP_df <- data.frame(EXP=sort(unique(x$EXP)))%>%
   mutate(values=as.numeric(gsub("B|M|R", "", EXP)))
EXP_df <- EXP_df%>%mutate(value2=ifelse(grepl("^M", EXP), values+12, values))

dd <- dd%>%mutate(Batch=str_sub(EXP, 1, 1))%>%
    left_join(EXP_df, by="EXP")%>%
    mutate(EXP=fct_reorder(EXP, value2))%>%as.data.frame()

##
### number of individual per library
dx <- x%>%distinct(sampleID, .keep_all=T)
dx <- dx%>%group_by(EXP)%>%summarise(nindi=n(),.groups="drop")%>%
   mutate(Batch=str_sub(EXP, 1, 1))%>%
   as.data.frame()
##
dx <- dx%>%left_join(EXP_df, by="EXP")%>%
    mutate(EXP=fct_reorder(EXP, value2))

p1 <- ggplot(dx)+
    geom_bar(stat="identity",  aes(x=EXP, y=nindi, fill=factor(Batch)))+
    ylab("#individual per library")+
    theme_bw()+
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=90, hjust=1, size=8))
##
png("./2_seurat.outs/Figure1.1_indi.png", width=550, height=420, res=120)
p1
dev.off()
          



p <- ggplot(dd)+
    geom_bar(stat="identity", aes(x=EXP, y=ncell, fill=factor(Batch)))+
    ggtitle("Number of cells")+
    theme_bw()+
    theme(legend.position="none",
       axis.title=element_blank(),
       axis.text.x=element_text(angle=90, hjust=1, size=8),
       plot.title=element_text(hjust=0.5))

png("./2_seurat.outs/Figure1.2_cells.png", width=550, height=420, res=120)
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

png("./2_seurat.outs/Figure1.3_reads.png", width=550, height=420, res=120)
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

png("./2_seurat.outs/Figure1.4_genes.png", width=550, height=420, res=120)
p
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
brainCV <- read.csv("BrainCV_cell.counts.csv")
brainCV <- brainCV%>%filter(USE==1)%>%
    mutate(Sample_ID2=toupper(Sample_ID2), comb=paste(Sample_ID2, Library, sep="_"))
###
brainCV2 <- brainCV%>%
   dplyr::select(EXP2=ExpectedPool, sampleID=Sample_ID2)%>%
   mutate(xvalue=as.numeric(gsub("^[BM]|R", "", EXP2)),
           xvalue=ifelse(grepl("^M", EXP2), xvalue+12, xvalue),
          sampleID2=paste(EXP2, sampleID, sep="_"))

brainCV2 <- brainCV%>%dplyr::select(sampleID=Sample_ID2, Library, Category)
dd2 <- dd%>%left_join(brainCV2, by="sampleID")
dd2 <- dd2%>%mutate(datasets=ifelse(grepl("^B", EXP), "Bannon", "Mash"))

dd2 <- dd <- dd%>%left_join(brainCV2, by="sampleID")
##
dd2 <- dd2%>%mutate(EXP=fct_reorder(EXP, yvalue), sampleID2=fct_reorder(sampleID2, xvalue))


### data for heatmap
dd3 <- dd2%>%group_by(sampleID)%>%
   mutate(Perc=ncell/sum(ncell))%>%as.data.frame()

### x labels
## sampleID <- dd2$sampleID
## names(sampleID) <- dd2$comb
p <- ggplot(dd3, aes(x=sampleID2, y=EXP, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_text(hjust=0, vjust=0.5, angle=-90, size=5),
         axis.text.y=element_text(size=8),
         axis.title=element_blank())

###
figfn <- "./2_seurat.outs/Figure2.1_heatmap.png"
png(figfn, width=1000, height=550, res=120)
print(p)
dev.off()


###
dd4 <- dd2%>%group_by(EXP)%>%
   mutate(Perc=ncell/sum(ncell))%>%as.data.frame()

### x labels
## sampleID <- dd2$sampleID
## names(sampleID) <- dd2$comb
p2 <- ggplot(dd4, aes(x=sampleID2, y=EXP, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_text(hjust=0, vjust=0.5, angle=-90, size=5),
         axis.text.y=element_text(size=8),
         axis.title=element_blank())

###
figfn <- "./2_seurat.outs/Figure2.2_heatmap.png"
png(figfn, width=1000, height=550, res=120)
print(p2)
dev.off()


###
### reads mapping to X and Y chromosome
##library(annotables)


### End


