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

## theme_set(theme_grey())

rm(list=ls())

###
### clustering analysis based on clean data

##
outdir <- "./3_cluster.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


sc <- read_rds("./2_seurat.outs/2_seurat.clean.rds")

## x <- sc@meta.data

###
sc <- sc%>%
   NormalizeData()%>%
   FindVariableFeatures(nfeatures=2000)%>%
   ScaleData()%>%
   RunPCA(npcs=100)

sc <- RunHarmony(sc, group.by.vars="EXP")
## sc <- sc%>%
##    RunUMAP(reduction="harmony", dims=1:50)%>%
##    FindNeighbors(dims=1:50, verbose=T)
sc <- sc%>%
   RunUMAP(reduction="harmony", dims=1:50)%>%
   FindNeighbors(reduction="harmony", dims=1:50, verbose=T)

###
sc2 <- sc%>%FindClusters(resolution=0.07, verbose=T)

opfn <- paste(outdir, "1_seurat.cluster14.rds", sep="")
write_rds(sc2, opfn)

###

all.genes <- rownames(sc)
sc <- sc%>%ScaleData(features=all.genes)
opfn <- "./3_cluster.outs/1_seurat.cluster14.rds"
write_rds(sc, opfn)
print(object.size(sc), units="Gb")




####
#### adjust resolution
sc <- read_rds("./3_cluster.outs/1_seurat.cluster14.rds")

mycol2 <- c("0"="#a6cee3", "1"="#1f78b4", "2"="#b2df8a", "3"="#33a02c",
   "4"="#fb9a99", "5"="#e31a1c", "6"="#fdbf6f", "7"="#ff7f00", "8"="#cab2d6",
    "9"="#6a3d9a","10"="#ffff99", "11"="#b15928") #, "12"="#8dd3c7", "13"="#8e0152", "14"="#d9d9d9")

## res_ls <- c(0.003, 0.01,  0.05, 0.1, 0.07, 0.2)
##
sc2 <- sc%>%FindClusters(resolution=0.1, verbose=T)
p0 <- DimPlot(sc2, reduction="umap", cols=mycol2, label=T, repel=T, pt.size=0.5)+
    theme_bw()
figfn <- "./3_cluster.outs/Figure0.6_cluster25.png"
png(figfn, width=520, height=480, res=120)
p0
dev.off()


###
x <- sc@meta.data%>%dplyr::select(seurat_clusters, sampleID)
cv <- read.csv("BrainCV_cell.counts.final.csv")%>%dplyr::filter(USE==1)%>%
   dplyr::select(sampleID, Library, Age, Race, Sex, BW_lbs, BW_g, pH,
                 CauseOfDeath, BrainRegion, Category, PC1, PC2, PC3)
x <- x%>%left_join(cv, by="sampleID")

x3 <- x%>%group_by(seurat_clusters)%>%summarise(ny=n(),.groups="drop")%>%mutate(percet=ny/sum(ny))

x2 <- x%>%group_by(seurat_clusters,Category)%>%summarise(ny=n(),.groups="drop")

x2 <- x2%>%group_by(seurat_clusters)%>%mutate(percet=ny/sum(ny))%>%as.data.frame()

x2%>%pivot_wider(id_cols=seurat_clusters, ,names_from=Category, values_from=percet)%>%as.data.frame()



#################
### UMAP plot ###
#################

sc <- read_rds("./3_cluster.outs/1_seurat.cluster.rds")

## ##
## mycol1 <- c("0"="#8dd3c7", "1"="#ffffb3", "2"="#bebada", "3"="#fb8072",
##   "4"="#80b1d3", "5"="#fdb462", "6"="#b3de69", "7"="#fccde5", "8"="#d9d9d9",
##   "9"="#bc80bd", "10"="#ccebc5", "11"="#ffed6f")

## mycol2 <- c("0"="#a6cee3", "1"="#1f78b4", "2"="#b2df8a", "3"="#33a02c",
##    "4"="#fb9a99", "5"="#e31a1c", "6"="#fdbf6f", "7"="#ff7f00", "8"="#cab2d6",
##     "9"="#6a3d9a","10"="#ffff99", "11"="#b15928", "12"="#8dd3c7", "13"="#8e0152", "14"="#d9d9d9") 

p0 <- DimPlot(sc, reduction="umap", cols=mycol2, label=T, repel=T, pt.size=0.5)+
    theme_bw()
figfn <- "./3_cluster.outs/Figure1.0_cluster.png"
png(figfn, width=520, height=480, res=120)
p0
dev.off()


## sc2 <- read_rds("./3_cluster.outs/1_seurat.cluster.rds")
umap <- as.data.frame(sc[["umap"]]@cell.embeddings)
x <- sc@meta.data
x <- x%>%mutate(Batch=ifelse(grepl("AKB", sampleID), "Bannon", "Mash"))
                            
 
df2 <- data.frame(UMAP_1=umap[,1],
   UMAP_2=umap[,2],
   seurat_clusters=x$seurat_clusters, Batch=x$Batch)
###
p <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(seurat_clusters)))+
   rasterise(geom_point(size=0.1),dpi=300)+
   scale_colour_manual(values=mycol2,
       guide=guide_legend(override.aes=list(size=2, ncol=1)))+
   ## guides(col=guide_legend(override.aes=list(size=2), ncol=1))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(1,"lines"))
###
figfn <- paste(outdir, "Figure1.1_cluster.png", sep="")
png(figfn, width=520, height=480, res=120)
p
dev.off()


### umap plot split by batch
p2 <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(seurat_clusters)))+
   rasterise(geom_point(size=0.1),dpi=300)+
       scale_colour_manual(values=mycol2,
          guide=guide_legend(override.aes=list(size=2, ncol=1)))+
   ## guides(col=guide_legend(override.aes=list(size=2), ncol=1))+
   facet_wrap(~factor(Batch), ncol=2)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(1,"lines"),
         strip.text=element_text(size=14))
###
figfn <- paste(outdir, "Figure1.2_cluster.png", sep="")
png(figfn, width=800, height=480, res=120)
p2
dev.off()


##################################################################
### find differentially expressed marker genes between cluster ###
##################################################################

sc <- read_rds("./3_cluster.outs/1_seurat.cluster14.rds")

##
## res=0.07, 14 cluster
brain.markers <- FindAllMarkers(sc, only.pos=T, min.pct=0.02, logfc.threshold=0.1)
opfn <- "./3_cluster.outs/2.0_cluster14_allmarker.rds"
write_rds(brain.markers, opfn)

##
x <- read_rds("./3_cluster.outs/2.0_cluster14_allmarker.rds")
x <- x%>%group_by(cluster)%>%arrange(desc(avg_log2FC), .by_group=T)%>%as.data.frame()
##
openxlsx::write.xlsx(x, "./3_cluster.outs/2.0_cluster14_allmarker.xlsx")x

##
brain.markers <- FindMarkers(sc, ident.1=7, only.pos=T, min.pct=0.02, logfc.threshold=0)
x <- brain.markers%>%arrange(desc(avg_log2FC))
openxlsx::write.xlsx(x, "./3_cluster.outs/2.0_cluster14_allmarker_7vsOthers.xlsx")


##
## x <- read_rds("./3_cluster.outs/2.0_cluster.marker.rds")
## x2 <- x%>%group_by(cluster)%>%slice_max(n=100, order_by=avg_log2FC)
## openxlsx::write.xlsx(x2, "./3_cluster.outs/2.0_top100_cluster.marker.xlsx")


###
### res=0.1, 17 clusters
sc2 <- sc%>%FindClusters(resolution=0.1, verbose=T)
brain.markers <- FindAllMarkers(sc2, only.pos=T, min.pct=0.2, logfc.threshold=0.1)
opfn <- "./3_cluster.outs/2.1_cluster17.marker.rds"
write_rds(brain.markers, opfn)
##
x <- read_rds("./3_cluster.outs/2.1_cluster17.marker.rds")
x <- x%>%group_by(cluster)%>%arrange(desc(avg_log2FC), .by_group=T)
##
openxlsx::write.xlsx(x, "./3_cluster.outs/2.1_cluster17.allmarker.xlsx")


###
### res=0.2, 25 cluster
sc2 <- sc%>%FindClusters(resolution=0.2, verbose=T)
brain.markers <- FindMarkers(sc2, ident.1=0, ident.2=1, min.pct=0.02, logfc.threshold=0)
opfn <- "./3_cluster.outs/2.2_cluster25.marker.rds"
write_rds(brain.markers, opfn)
###
x <- read_rds("./3_cluster.outs/2.2_cluster25.marker.rds")
x <- x%>%arrange(desc(abs(avg_log2FC)))
##
openxlsx::write.xlsx(x, "./3_cluster.outs/2.2_cluster25_0vs1.xlsx")




##################
### umap plots ###
##################

## geneList0 <- c("CTNNA3", "SLC24A2", "MBP", "MOBP", "PLP1", "MOG")
## geneList1 <- c("DOCK8", "C3", "CX3CR1", "LRRK1")
## geneList2 <- c("GFAP", "SLC4A4", "AQP4")
## geneList3 <- c("SOX6", "GRM7", "GRIK2", "DAB1","DCC", "GRM5")
## geneList4 <- c("RGS5", "SLC6A1")
## geneList5 <- c("SYT1", "SNAP25", "MEG3", "SLC18A2") #, "GRIA1", "GABBR2",
## geneList6 <- c("SKAP1", "THEMIS", "IL7R", "HLA-B") # "HLA-A")


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

outdir <- "./3_cluster.outs/"

### cluster 0
geneList0 <- c("CTNNA3", "MOG", "MOBP") 
fig_ls <- plotList(geneList0)
figfn <- paste(outdir, "Figure3.2.0_umap.cluster0.png", sep="")
png(figfn, width=720, height=300, res=120)
plot_grid(plotlist=fig_ls, ncol=3, align="hv")       
dev.off()


### cluster 1
geneList1 <- c("GFAP", "SLC4A4", "AQP4")
fig_ls <- plotList(geneList1)
figfn <- paste(outdir, "Figure3.2.1_umap.cluster1.png", sep="")
png(figfn, width=720, height=300, res=120)
plot_grid(plotlist=fig_ls,  ncol=3, align="hv")       
dev.off()


### cluster 2
geneList2 <- c("DOCK8", "C3", "CSF1R")
fig_ls <- plotList(geneList2)
figfn <- paste(outdir, "Figure3.2.2_umap.cluster2.png", sep="")
png(figfn, width=720, height=300, res=120)
plot_grid(plotlist=fig_ls, ncol=3, align="hv")       
dev.off()

 
### cluster 3 
geneList3 <- c("VCAN", "PDGFRA", "SOX6", "GRM7", "DCC", "GRM5")   
fig_ls <- plotList(geneList3)
figfn <- paste(outdir, "Figure3.2.3_umap.cluster3.png", sep="")
## png(figfn, width=1050, height=600, res=120)
png(figfn, width=720, height=600, res=120)
plot_grid(plotlist=fig_ls, ncol=3, align="hv")       
dev.off()

## geneList3 <- c("VCAN", "SOX6", "GRM7", "DCC", "GRIK2", "GRM5")
## fig_ls <- plotList(geneList3)
## figfn <- paste(outdir, "Figure3.2.3v2_umap.cluster3.png", sep="")
## png(figfn, width=720, height=600, res=120)
## plot_grid(plotlist=fig_ls, ncol=3, align="hv")       
## dev.off()

 
### cluster 4
geneList4 <- c("GALNTL6", "SYT1", "SYN3", "SYN2", "SNAP25",
    "GABRB2", "GABRB3", "GABBR2", "GAD1", "GAD2")
fig_ls <- plotList(geneList4)
figfn <- paste(outdir, "Figure3.2.4_umap.cluster4.png", sep="")
## png(figfn, width=900, height=600, res=120)
png(figfn, width=1050, height=600, res=120)
plot_grid(plotlist=fig_ls,  ncol=5, align="hv")       
dev.off()


### cluster 6
geneList6 <- c("RGS5", "CSPG4", "PECAM1", "KDR", "CDH5", "TEK")
fig_ls <- plotList(geneList6)
figfn <- paste(outdir, "Figure3.2.6_umap.cluster6.png", sep="")
png(figfn, width=720, height=600, res=120)
plot_grid(plotlist=fig_ls, ncol=3, align="hv")       
dev.off()



### cluster 7
geneList7 <- c("SLC18A2", "DDC", "TH", "SLC6A3",
    "NR4A2", "CALB1", "AGTR1", "ALDH1A1")
fig_ls <- plotList(geneList7)
figfn <- paste(outdir, "Figure3.2.7_umap.cluster7.png", sep="")
png(figfn, width=900, height=600, res=120)
plot_grid(plotlist=fig_ls, ncol=4, align="hv")       
dev.off()


### cluster 9
geneList6 <- c("SKAP1", "THEMIS", "IL7R")
fig_ls <- plotList(geneList6)
figfn <- paste(outdir, "Figure3.2.9_umap.cluster9.png", sep="")
png(figfn, width=720, height=300, res=120)
plot_grid(plotlist=fig_ls, ncol=3, align="hv")       
dev.off()





### specific marker gene expression
## geneList <- c("CTNNA3", "DOCK8", "GFAP", "GRM7", "RGS5", "SYT1", "IL7R")
## fig_ls <- plotList(geneList)
## figfn <- paste(outdir, "Figure2.2_umap.marker.png", sep="")
## png(figfn, width=850, height=400, res=120)
## plot_grid(plotlist=fig_ls,  ncol=4, align="hv")       
## dev.off()

## ###
## geneList <- c("GFAP", "OLR1", "GINS3", "TH",
##     "SLC6A3", "RGS5", "GAD2","CSF1R",
##     "MOG", "MOBP", "LGALS1")
## fig_ls <- plotList(geneList)
## figfn <- paste(outdir, "Figure2.3_umap.marker.png", sep="")
## png(figfn, width=850, height=700, res=120)
## plot_grid(plotlist=fig_ls, ncol=4, align="hv")       
## dev.off()

## ###
## geneList <- c("CTNNA3", "MOG", "MOBP", "DOCK8", "OLR1", "CSF1R", "GFAP", "GRM7", "SYT1", "TH", "SLC6A3")
## fig_ls <- plotList(geneList)
## figfn <- paste(outdir, "Figure2.4_umap.marker.png", sep="")
## png(figfn, width=850, height=700, res=120)
## plot_grid(plotlist=fig_ls, ncol=4, align="hv")       
## dev.off()

#################
### dot plots ###
#################
## ODC, astrocytes, microglia, DA, GABAergic, endothelial cells
geneList <- c("MOBP", "MBP", "PLP1", "CNP", 
   "GFAP", "AQP4", "SLC1A2",
   "CSF1R", "CX3CR1", "LRRK1", "DOCK8", 
   "VCAN", "PDGFRA",
   "GAD1", "GAD2", "GABRA2", "GABRA3", "GABRA4", "GABRA6",    
   "SLC18A2", "DDC", "TH", "SLC6A3", "NR4A2", "CALB1", "AGTR1",
   "RGS5", "CSPG4", "PECAM1", "KDR", "CDH5",
   "SKAP1", "THEMIS", "IL7R",
   "FOXJ1", "PIFO", "RABL2", "HSPH1")
 
##
MCls <- c("0"="ODC:0", "5"="ODC:5", "10"="ODC:10", "13"="ODC:13",
          "1"="astro:1", "2"="microglia:2", "11"="microglia:11",
          "3"="OPC:3",
          "4"="GABA:4", "7"="DaN:7",
          "6"="Endo:6", "8"="Endo:8",
          "9"="Tcell:9", "12"="Ependynal:12")
 
sc2 <- RenameIdents(sc, MCls)
 
##
p <- DotPlot(sc2, features=geneList, cols=c("blue", "red"), dot.scale=8)+
    theme(legend.title=element_text(size=10),
          axis.title.y=element_blank(),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(angle=45, hjust=1, size=8))

##
figfn <- "./3_cluster.outs/Figure4.1_dotplot.png"
png(figfn, width=950, height=500, res=120)
print(p)
dev.off()



################################################################
#### summary number of cells in clusters for each individual ###
################################################################

sc <- read_rds("./3_cluster.outs/1_seurat.cluster14.rds")

x <- sc@meta.data

x2 <- x%>%mutate(comb2=paste(EXP, sampleID, sep="_"))%>%
   group_by(comb2, seurat_clusters)%>%
   summarise(ncell=n(),.groups="drop")%>%as.data.frame()
x2 <- x2%>%mutate(sampleID=gsub(".*_", "", comb2),
    EXP=gsub("_.*", "", comb2))

x.DF <- x2%>%
    pivot_wider(id_cols=sampleID, names_from=seurat_clusters, names_prefix="cluster",
                values_from=ncell, values_fill=0)%>%
    as.data.frame()
x.DF$total <- rowSums(x.DF[,-1])

openxlsx::write.xlsx(x.DF, "./3_cluster.outs/3_summary_number_cells.xlsx", overwrite=T)

###
### barplots


mycol2 <- c("0"="#a6cee3", "1"="#1f78b4", "2"="#b2df8a", "3"="#33a02c",
   "4"="#fb9a99", "5"="#e31a1c", "6"="#fdbf6f", "7"="#ff7f00", "8"="#cab2d6",
    "9"="#6a3d9a","10"="#ffff99", "11"="#b15928", "12"="#8dd3c7", "13"="#8e0152") 

plotDF <- x2%>%group_by(sampleID)%>%mutate(prop=ncell/sum(ncell))
## %>%
##     mutate(EXP_value=as.numeric(gsub("M|B|R", "", EXP)),
##            EXP_value=ifelse(grepl("^M", EXP), as.numeric(EXP_value+12), as.numeric(EXP_value)),
##            comb2=fct_reorder(as.character(comb2), EXP_value))
###
p <- ggplot(plotDF, aes(x=as.character(comb2), y=prop, fill=seurat_clusters))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=mycol2)+
      ##guide=guide_legend(override.aes=list(size=0.5), ncol=2))+
    ylab("Proportion of #cells in cluster")+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.key.size=grid::unit(0.5, "cm"),
         axis.title.x=element_blank(),
         axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=5),
         plot.title=element_text(hjust=0.5))
 
figfn <- "./3_cluster.outs/Figure5.0_cells_barplot.png"
png(figfn, width=1000, height=420, res=120)
print(p)
dev.off()        
                                    



######################################################################
### summary the reads mapping to mitochondrion for each individual ###
######################################################################

sc <- read_rds("./3_cluster.outs/1_seurat.cluster14.rds")

meta <- sc@meta.data

###
### barplots
dd <- meta%>%mutate(comb2=paste(EXP, sampleID, sep="_"))%>%
     group_by(comb2)%>%summarise(ncell=n(), percent.mt=median(percent.mt),.groups="drop")

dd <- dd%>%mutate(sampleID=gsub(".*_", "", comb2),
                  EXP=gsub("_.*", "", comb2),
                  Batch=ifelse(grepl("AKB", sampleID), "Bannon", "Mash"))
##
p1 <- ggplot(dd, aes(x=comb2, y=percent.mt, fill=factor(EXP)))+
    geom_bar(stat="identity")+
    ggtitle("percentage of reads mapping to  mitochondrion")+
    theme_bw()+
    theme(legend.position="none",
          axis.title=element_blank(),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=5),
          plot.title=element_text(hjust=0.5))

figfn <- "./3_cluster.outs/Figure5.1_mt_barplot.png"
png(figfn, width=920, height=420, res=120)
print(p1)
dev.off()


###
### violin plot
meta2 <- meta%>%mutate(comb2=paste(EXP, sampleID, sep="_"))
p2 <- ggplot(meta2, aes(x=comb2, y=percent.mt, fill=factor(EXP)))+
   geom_violin(width=1)+
   ggtitle("percentage of reads mapping to  mitochondrion")+
   theme_bw()+
   theme(legend.position="none",
         axis.title=element_blank(),
         axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=5),
         plot.title=element_text(hjust=0.5))

figfn <- "./3_cluster.outs/Figure5.2_mt_violin.png"
png(figfn, width=920, height=420, res=120)
print(p2)
dev.off()    

## sc0 <- read_rds("./2_seurat.outs/1.1_seurat.SNG.rds")
                     
################
### heatmap ####
################
## sc <- read_rds("./3_cluster.outs/1_seurat.cluster14.rds")
## x <- read_rds("./3_cluster.outs/2.0_cluster.marker.rds")
## top10 <- x%>%group_by(cluster)%>%slice_max(n=10, order_by=avg_log2FC)%>%dplyr::pull(gene)%>%unique()


## mycol2 <- c("0"="#a6cee3", "1"="#1f78b4", "2"="#b2df8a", "3"="#33a02c",
##    "4"="#fb9a99", "5"="#e31a1c", "6"="#fdbf6f", "7"="#ff7f00", "8"="#cab2d6",
##     "9"="#6a3d9a","10"="#ffff99", "11"="#b15928")
## p <- DoHeatmap(sc, features=top10)+NoLegend()
 
## figfn <- "./3_cluster.outs/Figure4.1_feature.DoHeatmap.png"
## png(figfn, width=800, height=420, res=100)
## print(p)
## dev.off()



###############################
### marker fearures of umap ###
###############################

## sc <- read_rds("./3_cluster.outs/1_seurat.cluster.rds")


## ## gene <- rownames(sc)
## geneList <- c("SLC6A3", "NR4A2", "DRD2", "TH", "SLC18A2", "NEUROG2")
## fig_ls <- lapply(geneList, function(x0){
##    ### 
##    p <- VlnPlot(sc2, features=x0)&
##       theme_bw()+ 
##       theme(legend.position="none",
##           axis.title.y=element_text(size=8),
##           axis.title.x=element_blank(),
##           axis.text=element_text(size=8),
##           plot.title=element_text(hjust=0.5, size=8))
##    p
## })

## figfn <- paste(outdir, "Figure3.1_violin.png", sep="")
## png(figfn, width=520, height=620, res=120)
## plot_grid(plotlist=fig_ls, ncol=2, align="hv")       
## dev.off()




##########################################
### heatmap of clusters and cell types ###
##########################################

## scPred <- read_rds("/wsu/home/groups/bannonlab/sc_brains/Analysis/James/Merged_Reference_Bannon_Label_Transfer_06122022/BrainData_Labelled_06212022.rds")
## meta <- scPred@meta.data
## meta <- meta%>%dplyr::select(NEW_BARCODE, predicted.celltype)
## x <- sc@meta.data%>%dplyr::select(NEW_BARCODE, seurat_clusters)
## meta <- meta%>%left_join(x, by="NEW_BARCODE")


## df1 <- meta%>%
##    group_by(predicted.celltype, seurat_clusters)%>%
##    summarize(Freq=n(),.groups="drop")%>%
##    group_by(seurat_clusters)%>%
##    mutate(Perc=Freq/sum(Freq))

## p1 <- ggplot(df1, aes(x=seurat_clusters, y=predicted.celltype, fill=Perc))+
##    geom_tile()+
##    scale_fill_gradient("Fraction of cells",
##                        low="#ffffc8", high="#7d0025", na.value=NA)+
##    xlab("Cluster")+ylab("predicted.celltype")+
##    theme_bw()+theme(axis.text.x=element_text(hjust=0.5, size=10))

## ###
## png("./3_cluster.outs/Figure4.0_heatmap.png", width=700, height=600, res=120)
## print(p1)
## dev.off()
