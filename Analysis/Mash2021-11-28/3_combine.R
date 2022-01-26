##
library(Matrix)
library(tidyverse)
library(data.table)
library(Seurat)
library(parallel)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(harmony)
library(Rcpp)
##
library(ggrastr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

rm(list=ls())

outdir <- "./3_combine.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


############################
### combine two datasets ###
############################

### Bannon data
fn <- "../Bannon2021-11-28/2_seurat.outs/1.2_seurat.clean.rds"
sc <- read_rds(fn)

## x <- sc@meta.data
## dd <- x%>%group_by(SNG.BEST.GUESS,EXP)%>%summarise(ncell=n(), .groups="drop")%>%as.data.frame()

###
fn <- "./2_seurat.outs/1.3_seurat.clean2.rds"
sc2 <- read_rds(fn)

combine <- merge(sc, sc2, add.cell.ids=c("Bannon", "Mash"), project="sc-brain")
x <- combine@meta.data
x$NEW_BARCODE <- rownames(x)
x <- x%>%mutate(Batch=gsub("_.*", "", NEW_BARCODE))
###
combine <- AddMetaData(combine, x)
###
opfn <- paste(outdir, "1_seurat.combine.rds", sep="")
write_rds(combine, opfn)



###########################
### clustering analysis ###
###########################

fn <- paste(outdir, "1_seurat.combine.rds", sep="")
sc <- read_rds(fn)
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

opfn <- paste(outdir, "2.0_seurat.cluster.rds", sep="")
write_rds(sc2, opfn)


#################
### show umap ###
#################

sc2 <- read_rds("./3_combine.outs/2.0_seurat.cluster.rds")

umap <- as.data.frame(sc2[["umap"]]@cell.embeddings)
x <- sc2@meta.data
df2 <- data.frame(UMAP_1=umap[,1],
   UMAP_2=umap[,2],
   seurat_clusters=x$seurat_clusters, Batch=x$Batch)

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
   guides(col=guide_legend(override.aes=list(size=2), ncol=1))+
   facet_wrap(~factor(Batch), ncol=2)+
   theme_bw()+
   theme(## legend.position=c(0.8,0.15),
         legend.title=element_blank(),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(1,"lines"))
###
figfn <- paste(outdir, "Figure2.2_cluster.png", sep="")
png(figfn, width=650, height=380, res=120)
p2
dev.off()


###
### marker gene expression 
plotList <- function(geneList){
   ### 
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


###
###
sc2 <- read_rds("./3_combine.outs/2.0_seurat.cluster.rds")
###
geneList <- c("CTNNA3", "DOCK8", "GFAP", "GRM7", "RGS5", "SYT1", "IL7R")
fig_ls <- plotList(geneList)
figfn <- paste(outdir, "Figure2.3_umap.marker.png", sep="")
png(figfn, width=850, height=400, res=120)
plot_grid(fig_ls[[1]], fig_ls[[2]], fig_ls[[3]], fig_ls[[4]],
   fig_ls[[5]], fig_ls[[6]], fig_ls[[7]],  ncol=4, align="hv")       
dev.off()



########################################
### correct Bannon and Mash datasets ###
########################################


###################
### run harmony ###
###################

fn <- paste(outdir, "1_seurat.combine.rds", sep="")
sc <- read_rds(fn)
###
sc <- sc%>%
   NormalizeData()%>%
   FindVariableFeatures(nfeatures=2000)%>%
   ScaleData()%>%
   RunPCA(npcs=100) 
###
sc <- RunHarmony(sc, group.by.vars="Batch")

sc <- sc%>%RunUMAP(reduction="harmony", dims=1:50)%>%
   FindNeighbors(dims=1:50, verbose=T)
###
sc2 <- sc%>%FindClusters(resolution=0.02, verbose=T)

opfn <- paste(outdir, "3_seurat.cluster.rds", sep="")
write_rds(sc2, opfn)



#################
### show umap ###
#################

sc2 <- read_rds("./3_combine.outs/3_seurat.cluster.rds")

umap <- as.data.frame(sc2[["umap"]]@cell.embeddings)
x <- sc2@meta.data
df2 <- data.frame(UMAP_1=umap[,1],
   UMAP_2=umap[,2],
   seurat_clusters=x$seurat_clusters, Batch=x$Batch)

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
figfn <- paste(outdir, "Figure3.1_cluster.png", sep="")
png(figfn, width=520, height=480, res=120)
p1
dev.off()

### split by batch
p2 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(seurat_clusters)))+
   rasterise(geom_point(size=0.1),dpi=300)+
       ## scale_colour_manual(values=col0,
       ##     guide=guide_legend(override.aes=list(size=2)))+
   guides(col=guide_legend(override.aes=list(size=2), ncol=1))+
   facet_wrap(~factor(Batch), ncol=2)+
   theme_bw()+
   theme(## legend.position=c(0.8,0.15),
         legend.title=element_blank(),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(1,"lines"))
###
figfn <- paste(outdir, "Figure3.2_cluster.png", sep="")
png(figfn, width=650, height=380, res=120)
p2
dev.off()


###
###
### marker gene expression 
plotList <- function(geneList){
   ### 
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


###
###
sc2 <- read_rds("./3_combine.outs/3_seurat.cluster.rds")
###
geneList <- c("CTNNA3", "DOCK8", "GFAP", "GRM7", "RGS5", "SYT1", "IL7R")
fig_ls <- plotList(geneList)
figfn <- paste(outdir, "Figure3.3_umap.marker.png", sep="")
png(figfn, width=850, height=400, res=120)
plot_grid(fig_ls[[1]], fig_ls[[2]], fig_ls[[3]], fig_ls[[4]],
   fig_ls[[5]], fig_ls[[6]], fig_ls[[7]],  ncol=4, align="hv")       
dev.off()



######################################
### find differential marker genes ###
######################################

### rename cluster id 
sc2 <- read_rds("./3_combine.outs/3_seurat.cluster.rds")

CL <- Idents(sc2)
newCL <- rep(0, length(CL))
##
newCL[CL==2] <- 1
newCL[CL==1] <- 2
newCL[CL==3] <- 3
newCL[CL==4] <- 4
newCL[CL==5] <- 5
newCL[CL==6] <- 6
newCL[CL==7] <- 7
newCL[CL==8] <- 8
newCL[CL==9] <- 9

sc2$newCL <- newCL
##
opfn <- paste(outdir, "3.1_seurat.clusterNew.rds", sep="")
write_rds(sc2, opfn)


### remove cluster 8 and 9
sc <- read_rds("./3_combine.outs/3.1_seurat.clusterNew.rds")
x <- sc@meta.data
sc2 <- subset(sc, subset=newCL!=8&newCL!=9)


### find differentially expressed marker
brain.markers <- FindAllMarkers(sc2, only.pos=T, min.pct=0.25, logfc.threshold=0.25)

opfn <- paste(outdir, "3.2_cluster.marker.rds", sep="")
write_rds(brain.markers, opfn)

###
brain.markers <- read_rds("./3_combine.outs/3.2_cluster.marker.rds")

x <- brain.markers%>%group_by(cluster)%>%top_n(100, wt=avg_log2FC)

write_tsv(x, "./3_combine.outs/3.2_cluster.markers.top100.tsv")
