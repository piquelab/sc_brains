##
library(Matrix)
library(tidyverse)
library(data.table)
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

basefolder <- "/wsu/home/groups/bannonlab/sc_brains/counts_brain_2021Nov_Mash/"
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

fn <- paste(basefolder, "libList.txt", sep="")
libList <- read.table(fn)$V1

sc_ls <- mclapply(libList, function(ii){
    cat(ii, "\n")
    sc <- readCellRanger(basefolder, ii)
    sc
},mc.cores=14)

combined <- merge(sc_ls[[1]], sc_ls[-1], project="sc-brain")
x <- combined@meta.data
x$NEW_BARCODE <- rownames(x)

### add demuxlet results
meta <- read_rds("./1_demux2_output/1_demux_New.ALL.rds")
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


###############################
### summary of raw datasets ###
###############################

fn <- paste(outdir, "1.0_seurat.all.rds", sep="")
sc <- read_rds(fn)

x <- sc@meta.data

### droplet type
dd <- x%>%
   group_by(EXP, DROPLET.TYPE)%>%
   summarise(ncell=n(),.groups="drop")
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
   summarise(ncell=n(), nreads=mean(nCount_RNA), ngenes=mean(nFeature_RNA), .groups="drop")

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
p <- ggplot(dd)+
    geom_bar(stat="identity", aes(x=EXP, y=ngenes), fill="#1c9099")+
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
### heatmap for cells of sample in batch
dd <- x%>%group_by(SNG.BEST.GUESS, EXP)%>%summarise(ncell=n(), .groups="drop")
dd <- dd%>%mutate(sampleID=gsub("_.*", "", SNG.BEST.GUESS))
###
dd2 <- dd%>%group_by(SNG.BEST.GUESS)%>%slice(which.max(ncell))%>%
    mutate(comb=paste(EXP, sampleID, sep="_"))
dd2 <- dd2%>%as.data.frame()%>%dplyr::select(sampleID, comb)

dd <- dd%>%left_join(dd2, by="sampleID")


### data for heatmap
dd3 <- dd%>%group_by(SNG.BEST.GUESS)%>%
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
figfn <- "./2_seurat.outs/Figure0.5_heatmap.png"
png(figfn, width=800, height=320, res=100)
print(p)
dev.off()


###
dd4 <- dd%>%group_by(EXP)%>%
   mutate(Perc=ncell/sum(ncell))%>%as.data.frame()

### x labels
## sampleID <- dd2$sampleID
## names(sampleID) <- dd2$comb
p <- ggplot(dd4, aes(x=comb, y=EXP, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_blank(), #element_text(hjust=1, vjust=0.5, angle=90, size=6),
         axis.text.y=element_text(size=8),
         axis.title=element_blank())

###
figfn <- "./2_seurat.outs/Figure0.6_heatmap.png"
png(figfn, width=800, height=320, res=100)
print(p)
dev.off()


###################################
### summary of singlet datasets ###
###################################

### generate clean data
fn <- paste(outdir, "1.1_seurat.SNG.rds", sep="")
sc <- read_rds(fn)
x <- sc@meta.data%>%mutate(comb=paste(EXP, SNG.BEST.GUESS, sep="_"))

dd <- x%>%group_by(SNG.BEST.GUESS, EXP)%>%summarise(ncell=n(), .groups="drop")

###
dd2 <- dd%>%group_by(SNG.BEST.GUESS)%>%slice(which.max(ncell))%>%
    mutate(comb=paste(EXP, SNG.BEST.GUESS, sep="_"))
dd2 <- dd2%>%as.data.frame()%>%dplyr::select(SNG.BEST.GUESS, comb)

dd <- dd%>%left_join(dd2, by="sampleID")

##clean data
x2 <- x%>%
   dplyr::filter(EXP!="NYGC2-M13", comb%in%unique(dd2$comb))
dd3 <- x2%>%
   group_by(SNG.BEST.GUESS, EXP)%>%summarise(ncell=n(), .groups="drop")%>%as.data.frame()%>%
   dplyr::filter(ncell>20)
x2 <- x2%>%dplyr::filter(SNG.BEST.GUESS%in%unique(dd3$SNG.BEST.GUESS))

    
sc2 <- subset(sc, cells=rownames(x2))
opfn <- paste(outdir, "1.2_seurat.clean.rds", sep="")
write_rds(sc2, file=opfn)



### number of cells
sc2 <- read_rds("./2_seurat.outs/1.2_seurat.clean.rds")
x <- sc2@meta.data

### individuals across batch
dd <- x%>%
   group_by(SNG.BEST.GUESS, EXP)%>%
   summarise(ncell=n(), .groups="drop")
dd2 <- dd%>%group_by(EXP)%>%summarise(nind=n(),.groups="drop")

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
dd <- x%>%group_by(SNG.BEST.GUESS, EXP)%>%
   summarise(ncell=n(), nreads=mean(nCount_RNA), ngenes=mean(nFeature_RNA),.groups="drop")%>%
   as.data.frame() 
dd <- dd%>%mutate(sampleID=gsub("_.*", "", SNG.BEST.GUESS),
                  comb=paste(EXP,sampleID, sep="_"))
###

###
### number of cells per individual
p <- ggplot(dd)+
    geom_bar(stat="identity", aes(x=comb, y=ncell, fill=EXP))+
    ggtitle("Cells per individual")+
    theme_bw()+
    guides(fill=guide_legend(override.aes=list(size=0.5), ncol=1))+
    theme(legend.title=element_blank(),
       legend.text=element_text(size=8),   
       legend.key.size=grid::unit(0.5, "lines"),
       axis.title=element_blank(),
       axis.text.x=element_blank(),#element_text(angle=90, hjust=1, size=6),
       plot.title=element_text(hjust=0.5))

png("./2_seurat.outs/Figure1.1_cells.png", width=500, height=300, res=120)
p
dev.off()


###
### reads per individual
p <- ggplot(dd)+
    geom_bar(stat="identity", aes(x=comb, y=nreads, fill=EXP))+
    ggtitle("Reads per cell")+
    theme_bw()+
    theme(legend.position="none",
       axis.title=element_blank(),
       axis.text.x=element_blank(),#element_text(angle=90, hjust=1, size=8),
       plot.title=element_text(hjust=0.5))

png("./2_seurat.outs/Figure1.2_reads.png", width=420, height=320, res=120)
p
dev.off()


###
### number of genes per individual
p <- ggplot(dd)+
    geom_bar(stat="identity", aes(x=comb, y=ngenes, fill=EXP))+
    ggtitle("Genes per cell")+
    theme_bw()+
    theme(legend.position="none",
       axis.title=element_blank(),
       axis.text.x=element_blank(), #element_text(angle=90, hjust=1, size=8),
       plot.title=element_text(hjust=0.5))

png("./2_seurat.outs/Figure1.3_genes.png", width=420, height=320, res=120)
p
dev.off()


##################################################################
### heatmap show distribution of individual in different batch ###
##################################################################

fn <- paste(outdir, "1.1_seurat.SNG.rds", sep="")
sc <- read_rds(fn)
x <- sc@meta.data

dd <- x%>%group_by(SNG.BEST.GUESS, EXP)%>%summarise(ncell=n(), .groups="drop")
dd <- dd%>%mutate(sampleID=gsub("_.*", "", SNG.BEST.GUESS))
###
dd2 <- dd%>%group_by(SNG.BEST.GUESS)%>%slice(which.max(ncell))%>%
    mutate(comb=paste(EXP, sampleID, sep="_"))
dd2 <- dd2%>%as.data.frame()%>%dplyr::select(sampleID, comb)

dd <- dd%>%left_join(dd2, by="sampleID")


### data for heatmap
dd3 <- dd%>%group_by(SNG.BEST.GUESS)%>%
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
figfn <- "./2_seurat.outs/Figure1.5_heatmap.png"
png(figfn, width=800, height=320, res=100)
print(p)
dev.off()


###
dd4 <- dd%>%group_by(EXP)%>%
   mutate(Perc=ncell/sum(ncell))%>%as.data.frame()

### x labels
## sampleID <- dd2$sampleID
## names(sampleID) <- dd2$comb
p <- ggplot(dd4, aes(x=comb, y=EXP, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_blank(), #element_text(hjust=1, vjust=0.5, angle=90, size=6),
         axis.text.y=element_text(size=8),
         axis.title=element_blank())

###
figfn <- "./2_seurat.outs/Figure1.6_heatmap.png"
png(figfn, width=800, height=320, res=100)
print(p)
dev.off()


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




###########################
### clustering analysis ###
###########################
 
sc <- read_rds("./2_seurat.outs/1.2_seurat.clean.rds")
x <- sc@meta.data
## dd <- x%>%
##    group_by(SNG.BEST.GUESS, EXP)%>%summarise(ncell=n(), .groups="drop")%>%as.data.frame()
##    ## dplyr::filter(ncell>20)


###
sc <- sc%>%
   NormalizeData()%>%
   FindVariableFeatures(nfeatures=2000)%>%
   ScaleData()
###
sc <- sc%>%
   RunPCA(npcs=100)%>%
   RunUMAP(dims=1:50)%>%
   FindNeighbors(dims=1:50, verbose=T)
###
sc2 <- sc%>%FindClusters(resolution=0.02, verbose=T)

opfn <- paste(outdir, "2.1_seurat.cluster.rds", sep="")
write_rds(sc2, opfn)



#################
### UMAP plot ###
#################

sc2 <- read_rds("./2_seurat.outs/2.1_seurat.cluster.rds")

CL <- Idents(sc2)

newCL <- rep(0, length(CL))
newCL[CL==3] <- 1
newCL[CL==1] <- 2
newCL[CL==4] <- 3
newCL[CL==6] <- 4
newCL[CL==2|CL==5] <- 5
newCL[CL==7] <- 6

sc2$newCL <- newCL
### output
opfn <- paste(outdir, "2.1_seurat.cluster.rds", sep="")
write_rds(sc2, opfn)

###
###
## Idents(sc2) <- as.numeric(sc2$newCL)
## ## col0 <- c("0"=, "1"=, "2"=, "3"=, "4"=, "5"=, "6"="")
## p0 <- DimPlot(object=sc2, label=T, raster=T)+
##    guides(col=guide_legend(override.aes=list(size=2), ncol=1)) &
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(1,"lines"))
## ### 
## figfn <- paste(outdir, "Figure2.0_cluster.png", sep="")
## png(figfn, width=520, height=480, res=120)
## p0
## dev.off()


umap <- as.data.frame(sc2[["umap"]]@cell.embeddings)
x <- sc2@meta.data
df2 <- data.frame(UMAP_1=umap[,1],
   UMAP_2=umap[,2],
   seurat_clusters=x$newCL, Batch=x$EXP)

### UMAP colored by cluster
p1 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(seurat_clusters)))+
   rasterise(geom_point(size=0.1),dpi=300)+
       ## scale_colour_manual(values=col0,
       ##     guide=guide_legend(override.aes=list(size=2)))+
   guides(col=guide_legend(override.aes=list(size=2), ncol=1))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(1,"lines"))
###
figfn <- paste(outdir, "Figure2.1_cluster.png", sep="")
png(figfn, width=520, height=480, res=120)
p1
dev.off()

### split by batch
p2 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(seurat_clusters)))+
   rasterise(geom_point(size=0.1),dpi=300)+
       ## scale_colour_manual(values=col0,
       ##     guide=guide_legend(override.aes=list(size=2)))+
   guides(col=guide_legend(override.aes=list(size=2), ncol=2))+
   facet_wrap(~factor(Batch), ncol=5)+
   theme_bw()+
   theme(legend.position=c(0.8,0.15),
         legend.title=element_blank(),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(1,"lines"))
###
figfn <- paste(outdir, "Figure2.2_cluster.png", sep="")
png(figfn, width=720, height=500, res=120)
p2
dev.off()
    

### further clean data
sc <- read_rds("./2_seurat.outs/1.2_seurat.clean.rds")
x <- sc@meta.data
dd <- x%>%
   group_by(SNG.BEST.GUESS, EXP)%>%summarise(ncell=n(), .groups="drop")%>%as.data.frame()
d2 <- dd%>%dplyr::filter(ncell>500)

x2 <- x%>%dplyr::filter(SNG.BEST.GUESS%in%unique(d2$SNG.BEST.GUESS))

sc2 <- subset(sc, cells=rownames(x2))
opfn <- paste(outdir, "1.3_seurat.clean2.rds", sep="")
write_rds(sc2, opfn)




###############################
### marker gene expression  ###
###############################

sc <- read_rds("./2_seurat.outs/2.1_seurat.cluster.rds")

##################
### umap plots ###
##################

geneList0 <- c("CTNNA3", "SLC24A2", "MBP", "MOBP", "PLP1", "MOG")
geneList1 <- c("DOCK8", "C3", "CX3CR1", "LRRK1")
geneList2 <- c("GFAP", "SLC4A4", "AQP4")
geneList3 <- c("SOX6", "GRM7", "GRIK2", "DAB1","DCC", "GRM5")
geneList4 <- c("RGS5", "SLC6A1")
geneList5 <- c("SYT1", "SNAP25", "MEG3", "SLC18A2") #, "GRIA1", "GABBR2",
geneList6 <- c("SKAP1", "THEMIS", "IL7R", "HLA-B") # "HLA-A")


## plotList function
plotList <- function(geneList){
    
   fig_ls <- lapply(geneList, function(x0){
   ###
      p <- FeaturePlot(sc, features=x0)+
         scale_color_gradient(low="lightgrey", high="blue")+
         ggtitle(bquote(~italic(.(x0))))&
         theme_bw()+ 
         theme(legend.position="none",
           axis.title=element_text(size=8),
           axis.text=element_text(size=8),
           plot.title=element_text(hjust=0.5, size=8))
      p
   })    
   fig_ls
}    


### cluster 0
fig_ls <- plotList(geneList0)
figfn <- paste(outdir, "Figure3.2.0_umap.cluster0.png", sep="")
png(figfn, width=720, height=400, res=120)
plot_grid(fig_ls[[1]], fig_ls[[2]], fig_ls[[3]], fig_ls[[4]],
   fig_ls[[5]], fig_ls[[6]], ncol=3, align="hv")       
dev.off()

### cluster 1
fig_ls <- plotList(geneList1)
figfn <- paste(outdir, "Figure3.2.1_umap.cluster1.png", sep="")
png(figfn, width=480, height=400, res=120)
plot_grid(fig_ls[[1]], fig_ls[[2]], fig_ls[[3]], fig_ls[[4]], ncol=2, align="hv")       
dev.off()

### cluster 2
fig_ls <- plotList(geneList2)
figfn <- paste(outdir, "Figure3.2.2_umap.cluster2.png", sep="")
png(figfn, width=600, height=280, res=120)
plot_grid(fig_ls[[1]], fig_ls[[2]], fig_ls[[3]],  ncol=3, align="hv")       
dev.off()

### cluster 3
fig_ls <- plotList(geneList3)
figfn <- paste(outdir, "Figure3.2.3_umap.cluster3.png", sep="")
png(figfn, width=720, height=400, res=120)
plot_grid(fig_ls[[1]], fig_ls[[2]], fig_ls[[3]], fig_ls[[4]],
   fig_ls[[5]], fig_ls[[6]], ncol=3, align="hv")       
dev.off()


### cluster 4
fig_ls <- plotList(geneList4)
figfn <- paste(outdir, "Figure3.2.4_umap.cluster4.png", sep="")
png(figfn, width=420, height=280, res=120)
plot_grid(fig_ls[[1]], fig_ls[[2]],  ncol=2, align="hv")       
dev.off()


### cluster 5
fig_ls <- plotList(geneList5)
figfn <- paste(outdir, "Figure3.2.5_umap.cluster5.png", sep="")
png(figfn, width=480, height=400, res=120)
plot_grid(fig_ls[[1]], fig_ls[[2]], fig_ls[[3]], fig_ls[[4]], ncol=2, align="hv")       
dev.off()


### cluster 6
fig_ls <- plotList(geneList6)
figfn <- paste(outdir, "Figure3.2.6_umap.cluster6.png", sep="")
png(figfn, width=480, height=400, res=120)
plot_grid(fig_ls[[1]], fig_ls[[2]], fig_ls[[3]], fig_ls[[4]], ncol=2, align="hv")       
dev.off()



### specific marker gene expression
geneList <- c("CTNNA3", "DOCK8", "GFAP", "GRM7", "RGS5", "SYT1", "IL7R")
fig_ls <- plotList(geneList)
figfn <- paste(outdir, "Figure3.2_umap.marker.png", sep="")
png(figfn, width=850, height=400, res=120)
plot_grid(fig_ls[[1]], fig_ls[[2]], fig_ls[[3]], fig_ls[[4]],
   fig_ls[[5]], fig_ls[[6]], fig_ls[[7]],  ncol=4, align="hv")       
dev.off()


