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
outdir <- "./3_cluster.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


sc <- read_rds("./2_seurat.outs/2.2_seurat.clean.rds")
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

opfn <- paste(outdir, "1.1_seurat.cluster.rds", sep="")
write_rds(sc2, opfn)


### UMAP plot
sc2 <- read_rds("./2_seurat.outs/1.1_seurat.cluster.rds")
umap <- as.data.frame(sc2[["umap"]]@cell.embeddings)
x <- sc2@meta.data
x <- x%>%mutate(EXP2=gsub("NYGC.-", "", EXP),
                Batch=ifelse(grepl("B", EXP2), "Bannon", "Mash"))
                            

df2 <- data.frame(UMAP_1=umap[,1],
   UMAP_2=umap[,2],
   seurat_clusters=x$seurat_clusters, Batch=x$Batch)
###
p <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(seurat_clusters)))+
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
p
dev.off()


### umap plot split by batch
p <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(seurat_clusters)))+
   rasterise(geom_point(size=0.1),dpi=300)+
       ## scale_colour_manual(values=col0,
       ##     guide=guide_legend(override.aes=list(size=2)))+
   ## guides(col=guide_legend(override.aes=list(size=2), ncol=1))+
   facet_wrap(~factor(Batch), ncol=2)+
   theme_bw()+
   theme(legend.position="none")## legend.title=element_blank(),
         ## legend.background=element_blank(),
         ## legend.box.background=element_blank(),
         ## legend.key.size=grid::unit(1,"lines"))
###
figfn <- paste(outdir, "Figure2.2_cluster.png", sep="")
png(figfn, width=700, height=480, res=120)
p
dev.off()



#################################################
### Finding differentially expressed features ###
#################################################

sc <- read_rds("./3_seurat.outs/1.1_seurat.cluster.rds") 
brain.markers <- FindAllMarkers(sc, only.pos=T, min.pct=0.25, logfc.threshold=0.25)

opfn <- paste(outdir, "2.0_cluster.marker.rds", sep="")
write_rds(brain.markers, opfn)

###
brain.markers <- read_rds("./2_seurat.outs/2.0_cluster.marker.rds")

x <- brain.markers%>%group_by(cluster)%>%top_n(100, wt=avg_log2FC)

write.csv(x, "./2_seurat.outs/2.0_cluster.markers.top100.csv", row.names=F)


#######################
### marker fearures ###
#######################

sc <- read_rds("./2_seurat.outs/2.1_seurat.cluster.rds")

## gene <- rownames(sc)
geneList <- c("SLC6A3", "NR4A2", "DRD2", "TH", "SLC18A2", "NEUROG2")
fig_ls <- lapply(geneList, function(x0){
   ### 
   p <- VlnPlot(sc, features=x0)&
      theme_bw()+ 
      theme(legend.position="none",
          axis.title.y=element_text(size=8),
          axis.title.x=element_blank(),
          axis.text=element_text(size=8),
          plot.title=element_text(hjust=0.5, size=8))
   p
})

figfn <- paste(outdir, "Figure3.1_violin.png", sep="")
png(figfn, width=520, height=620, res=120)
plot_grid(fig_ls[[1]], fig_ls[[2]], fig_ls[[3]], fig_ls[[4]],
   fig_ls[[5]], fig_ls[[6]], ncol=2, align="hv")       
dev.off()


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
