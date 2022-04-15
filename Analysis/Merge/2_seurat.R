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
theme_set(theme_grey())

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
},mc.cores=7)

combined <- merge(sc_ls[[1]], sc_ls[-1], project="sc-brain")
x <- combined@meta.data
x$NEW_BARCODE <- rownames(x)

### add demuxlet results
#### 13,068 cells not shown in the demuxlet, B3R, B4R and M13  
meta <- read_rds("./1_demux2_v2_output/1_demux_New.ALL.rds")
meta2 <- meta%>%dplyr::select("DROPLET.TYPE", "SNG.BEST.GUESS", "NEW_BARCODE", "EXP")
x <- x%>%left_join(meta2, by="NEW_BARCODE")
rownames(x) <- x$NEW_BARCODE

combined <- AddMetaData(combined, x)

opfn <- paste(outdir, "1.0_seurat.all.rds", sep="")
write_rds(combined, opfn)


### clean data
fn <- paste(outdir, "1.0_seurat.all.rds", sep="")
sc <- read_rds(fn)
sc2 <- subset(sc, subset=DROPLET.TYPE=="SNG")
###
opfn <- paste(outdir, "1.1_seurat.SNG.rds", sep="")
write_rds(sc2, opfn)




################################
### merge different datasets ###
################################

### 
sc1 <- read_rds("./2_seurat.outs/1.0_seurat.all.rds")
x <- sc1@meta.data%>%drop_na(EXP)
sc1 <- subset(sc1, cells=x$NEW_BARCODE)
    
###
prefix <- "../Bannon_2022Mar/"
fn <- paste(prefix, "2_seurat.outs/1.0_seurat.all.rds", sep="")
sc2 <- read_rds(fn)
expSel <- c("NYGC3-B7", "NYGC3-B8", "NYGC3-B9", "NYGC3-B10", "NYGC3-B11", "NYGC3-B12")
sc2 <- subset(sc2, subset=EXP%in%expSel)

### M1-M14
prefix <- "../Mash2021-11-28/"
fn <- paste(prefix, "2_seurat.outs/1.0_seurat.all.rds", sep="")
sc3 <- read_rds(fn)
sc3 <- subset(sc3, subset=EXP!="NYGC2-M13")

### M15-M20
prefix <- "../Mash_2022Mar/"
fn <- paste(prefix, "2_seurat.outs/1.0_seurat.all.rds", sep="")
sc4 <- read_rds(fn)
sc4 <- subset(sc4, subset=EXP!="NYGC3-M13")

###
combined <- merge(x=sc1, y=list(sc2, sc3, sc4))
opfn <- paste(outdir, "2.0_seurat.all.rds", sep="")
write_rds(combined, opfn)

### 
combined2 <- subset(combined, subset=DROPLET.TYPE=="SNG")
opfn <- paste(outdir, "2.1_seurat.SNG.rds", sep="")
write_rds(combined2, opfn)


    
###############################
### summary of all datasets ###
###############################

fn <- paste(outdir, "2.1_seurat.SNG.rds", sep="")
sc <- read_rds(fn)

x <- sc@meta.data

### droplet type
dd <- x%>%
   group_by(EXP, DROPLET.TYPE)%>%
   summarise(ncell=n(),.groups="drop")

dd <- dd%>%drop_na(EXP) ## I will check NA values

d2 <- dd%>%group_by(EXP)%>%mutate(perc=ncell/sum(ncell))%>%ungroup()
d2%>%filter(DROPLET.TYPE=="DBL")


p <- ggplot(dd)+
   geom_bar(stat="identity", position=position_fill(reverse=F),
        aes(x=EXP, y=ncell, fill=DROPLET.TYPE))+
   theme_bw()+
   theme(legend.title=element_blank(),
         axis.title=element_blank(),
         axis.text.x=element_text(angle=90, hjust=1, size=8))

png("./2_seurat.outs/Figure0.1_droplet.png", width=420, height=420, res=120)
p
dev.off()


### number of cells
dd <- x%>%
   group_by(EXP)%>%
   summarise(ncell=n(), .groups="drop")
dd <- dd%>%drop_na(EXP)
p <- ggplot(dd)+
    geom_bar(stat="identity", aes(x=EXP, y=ncell), fill="#1c9099")+
    ggtitle("Number of cells")+
    theme_bw()+
    theme(legend.title=element_blank(),
       axis.title=element_blank(),
       axis.text.x=element_text(angle=90, hjust=1, size=8),
       plot.title=element_text(hjust=0.5))

png("./2_seurat.outs/Figure0.2_cells.png", width=400, height=500, res=120)
p
dev.off()


### reads per cell
dd <- x%>%
   group_by(EXP)%>%
   summarise(nreads=mean(nCount_RNA), .groups="drop")
dd <- dd%>%drop_na(EXP)
p <- ggplot(dd)+
    geom_bar(stat="identity", aes(x=EXP, y=nreads), fill="#1c9099")+
    ggtitle("Reads per cell")+
    theme_bw()+
    theme(legend.title=element_blank(),
       axis.title=element_blank(),
       axis.text.x=element_text(angle=90, hjust=1, size=8),
       plot.title=element_text(hjust=0.5))

png("./2_seurat.outs/Figure0.3_reads.png", width=400, height=500, res=120)
p
dev.off()


###
### reads per cell
dd <- x%>%
   group_by(EXP)%>%
   summarise(nFeature=mean(nFeature_RNA), .groups="drop")
dd <- dd%>%drop_na(EXP)
p <- ggplot(dd)+
    geom_bar(stat="identity", aes(x=EXP, y=nFeature), fill="#1c9099")+
    ggtitle("Genes per cell")+
    theme_bw()+
    theme(legend.title=element_blank(),
       axis.title=element_blank(),
       axis.text.x=element_text(angle=90, hjust=1, size=8),
       plot.title=element_text(hjust=0.5))

png("./2_seurat.outs/Figure0.4_genes.png", width=400, height=500, res=120)
p
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

x2 <- x%>%select(sampleID=SNG.BEST.GUESS, EXP)
dd <- x2%>%group_by(sampleID, EXP)%>%
   summarise(ncell=n(), .groups="drop")


###
dd2 <- dd%>%group_by(sampleID)%>%top_n(1, wt=ncell)%>%
    mutate(comb=paste(EXP, sampleID, sep="_"))
dd2 <- dd2%>%as.data.frame()%>%dplyr::select(sampleID, comb)

dd <- dd%>%left_join(dd2, by="sampleID")
opfn <- "./2_seurat.outs/sample_EXP.comb.rds"
write_rds(dd, opfn)

 
## dd <- dd%>%drop_na(EXP)

### data for heatmap
dd3 <- dd%>%group_by(sampleID)%>%
   mutate(Perc=ncell/sum(ncell))%>%as.data.frame()

### x labels
## sampleID <- dd2$sampleID
## names(sampleID) <- dd2$comb
p <- ggplot(dd3, aes(x=comb, y=EXP, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_blank(), #element_text(hjust=1, vjust=0.5, angle=90, size=6),
         axis.text.y=element_text(size=8),
         axis.title=element_blank())

###
figfn <- paste(outdir, "Figure2.5_heatmap.png", sep="")
png(figfn, width=600, height=300, res=100)
print(p)
dev.off()



##################
### clean data ###
##################

sc <- read_rds("./2_seurat.outs/2.1_seurat.SNG.rds")
x <- sc@meta.data%>%mutate(comb=paste(EXP, "_", SNG.BEST.GUESS, sep=""))
dd <- read_rds("./2_seurat.outs/sample_EXP.comb.rds")
d2 <- dd%>%filter(ncell>100)
combSel <- unique(d2$comb)
cellSel <- x%>%filter(comb%in%combSel)%>%dplyr::pull(NEW_BARCODE)
###
sc2 <- subset(sc, cells=cellSel)
opfn <- "./2_seurat.outs/2.2_seurat.clean.rds"
write_rds(sc2, opfn)


x <- sc2@meta.data


###################################
### summary of singlet datasets ###
###################################


fn <- paste(outdir, "1.1_seurat.SNG.rds", sep="")
sc <- read_rds(fn)
x <- sc@meta.data%>%mutate(comb=paste(SNG.BEST.GUESS, EXP, sep="_"))
x2 <- x%>%dplyr::filter(!SNG.BEST.GUESS%in%c("AKB23", "AKB24", "AKB39", "AKB49", "AKB56"),
   EXP!="NYGC2-B4R",
   !comb%in%c("AKB14_NYGC2-B1R", "AKB20_NYGC2-B6R", "AKB30_NYGC2-B2R", "AKB32_NYGC2-B6R",
             "AKB32_NYGC2-B1R", "AKB35_NYGC2-B1R"))

sc2 <- subset(sc, cells=rownames(x2))
opfn <- paste(outdir, "1.2_seurat.clean.rds", sep="")
write_rds(sc2, file=opfn)


### number of cells

sc2 <- read_rds("./2_seurat.outs/1.2_seurat.clean.rds")
x <- sc2@meta.data

dd <- x%>%
   group_by(SNG.BEST.GUESS, EXP)%>%
   summarise(ncell=n(), .groups="drop")
dd2 <- dd%>%group_by(EXP)%>%summarise(nind=n(),.groups="drop")


### individuals across batch
p <- ggplot(dd2)+
    geom_bar(stat="identity", aes(x=EXP, y=nind), fill="#1c9099")+
    ggtitle("Number of individual")+
    theme_bw()+
    theme(legend.title=element_blank(),
       axis.title=element_blank(),
       axis.text.x=element_text(angle=90, hjust=1, size=8),
       plot.title=element_text(hjust=0.5))

png("./2_seurat.outs/Figure1.0_ind.png", width=450, height=420, res=120)
p
dev.off()


###
###
dd <- x%>%group_by(SNG.BEST.GUESS)%>%
   summarise(ncell=n(), nreads=mean(nCount_RNA), ngenes=mean(nFeature_RNA),.groups="drop")
d2 <- dd%>%filter(ncell>200)
###

###
### number of cells per individual
p <- ggplot(dd)+
    geom_bar(stat="identity", aes(x=SNG.BEST.GUESS, y=ncell, fill=EXP))+
    ggtitle("Cells per individual")+
    theme_bw()+
    theme(legend.title=element_blank(),
       axis.title=element_blank(),
       axis.text.x=element_text(angle=90, hjust=1, size=8),
       plot.title=element_text(hjust=0.5))

png("./2_seurat.outs/Figure1.1_cells.png", width=600, height=400, res=120)
p
dev.off()


###
### reads per individual
p <- ggplot(dd)+
    geom_bar(stat="identity", aes(x=SNG.BEST.GUESS, y=nreads, fill=EXP))+
    ggtitle("Reads per cell")+
    theme_bw()+
    theme(legend.position="none",
       axis.title=element_blank(),
       axis.text.x=element_text(angle=90, hjust=1, size=8),
       plot.title=element_text(hjust=0.5))

png("./2_seurat.outs/Figure1.2_reads.png", width=480, height=420, res=120)
p
dev.off()


###
### number of genes per individual
p <- ggplot(dd)+
    geom_bar(stat="identity", aes(x=SNG.BEST.GUESS, y=ngenes, fill=EXP))+
    ggtitle("Genes per cell")+
    theme_bw()+
    theme(legend.position="none",
       axis.title=element_blank(),
       axis.text.x=element_text(angle=90, hjust=1, size=8),
       plot.title=element_text(hjust=0.5))

png("./2_seurat.outs/Figure1.3_genes.png", width=480, height=420, res=120)
p
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
## figfn <- "./2_seurat.outs/Figure1.4.0_heatmap.png"
## png(figfn, width=480, height=280, res=100)
## print(p)
## dev.off()


###
### clean data

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
## figfn <- "./2_seurat.outs/Figure1.4.1_heatmap.png"
## png(figfn, width=480, height=280, res=100)
## print(p)
## dev.off()


####################################################################
### heatmap show distribution of individual in different batches ###
####################################################################

fn <- paste(outdir, "1.1_seurat.SNG.rds", sep="")
sc <- read_rds(fn)
x <- sc@meta.data

dd <- x%>%group_by(SNG.BEST.GUESS, EXP)%>%summarise(ncell=n(), .groups="drop")
## dd <- dd%>%mutate(sampleID=gsub("_.*", "", SNG.BEST.GUESS))
## ###
## dd2 <- dd%>%group_by(SNG.BEST.GUESS)%>%slice(which.max(ncell))%>%
##     mutate(comb=paste(EXP, sampleID, sep="_"))
## dd2 <- dd2%>%as.data.frame()%>%dplyr::select(sampleID, comb)

## dd <- dd%>%left_join(dd2, by="sampleID")


### data for heatmap
dd3 <- dd%>%group_by(SNG.BEST.GUESS)%>%
   mutate(Perc=ncell/sum(ncell))%>%as.data.frame()

### x labels
## sampleID <- dd2$sampleID
## names(sampleID) <- dd2$comb
p <- ggplot(dd3, aes(x=SNG.BEST.GUESS, y=EXP, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_text(hjust=1, vjust=0.5, angle=90, size=6),
         axis.text.y=element_text(size=8),
         axis.title=element_blank())

###
figfn <- "./2_seurat.outs/Figure1.5_heatmap.png"
png(figfn, width=600, height=280, res=100)
print(p)
dev.off()


###
dd4 <- dd%>%group_by(EXP)%>%
   mutate(Perc=ncell/sum(ncell))%>%as.data.frame()

### x labels
## sampleID <- dd2$sampleID
## names(sampleID) <- dd2$comb
p <- ggplot(dd4, aes(x=SNG.BEST.GUESS, y=EXP, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_text(hjust=1, vjust=0.5, angle=90, size=6),
         axis.text.y=element_text(size=8),
         axis.title=element_blank())

###
figfn <- "./2_seurat.outs/Figure1.6_heatmap.png"
png(figfn, width=600, height=280, res=100)
print(p)
dev.off()





###########################
### clustering analysis ###
###########################
 

