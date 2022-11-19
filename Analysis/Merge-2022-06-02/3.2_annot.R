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
outdir <- "./3.2_annot.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

## query data
sc <- read_rds("./2_seurat.outs/2_seurat.clean.rds")
sc <- sc%>%
   NormalizeData()%>%
   FindVariableFeatures(nfeatures=2000)%>%
   ScaleData()%>%
   RunPCA(npcs=100)


############################
### First reference data ###
############################

## load("/wsu/home/groups/bannonlab/sc_brains/Analysis/ReferenceGSE140231/Part3_Label_Clusters_06072022/Harmony_Clustered_UMAP20_Res04_substantia_nigra__big_06122022.RData")
## ref <- substantia_nigra_big2
## write_rds(ref, "refdata1_nc2020.rds")
ref <- read_rds("refdata1_nc2020.rds")
ref <- ref%>%
   NormalizeData()%>%
   FindVariableFeatures(nfeatures=2000)%>%
   ScaleData()%>%
   RunPCA(npcs=100)


### transfer cell-type
anchors <- FindTransferAnchors(reference=ref,  query=sc, dims=1:50, reference.reduction="pca")
pred <- TransferData(anchorset=anchors,  refdata=ref$CellType, dims=1:50)
pred$NEW_BARCODE <- rownames(pred)

opfn <- "./3.2_annot.outs/1_pred.rds"
write_rds(pred, opfn)



#########################
### umap of cell-type ###
#########################

sc <- read_rds("./3_cluster.outs/1_seurat.cluster14.rds")
x <- sc@meta.data

pred <- read_rds("./3.2_annot.outs/1_pred.rds")
pred2 <- pred%>%dplyr::select(NEW_BARCODE, predicted.id)

x2 <- x%>%left_join(pred2, by="NEW_BARCODE")
rownames(x2) <- x2$NEW_BARCODE

sc2 <- AddMetaData(sc, metadata=x2)
Idents(sc2) <- sc2$predicted.id

## ##
## mycol1 <- c("0"="#8dd3c7", "1"="#ffffb3", "2"="#bebada", "3"="#fb8072",
##   "4"="#80b1d3", "5"="#fdb462", "6"="#b3de69", "7"="#fccde5", "8"="#d9d9d9",
##   "9"="#bc80bd", "10"="#ccebc5", "11"="#ffed6f")

## mycol2 <- c("0:ODC"="#a6cee3", "1:ODC"="#1f78b4", "2:ODC"="#b2df8a", "3:Astrocyte"="#33a02c",
##    "4:Microglia"="#fb9a99", "5:OPC"="#e31a1c", "6:Astrocyte"="#fdbf6f", "7:Dopamine Neurons"="#ff7f00",
##    "8:Endothelial"="#cab2d6", "9:GABA"="#6a3d9a")
##"10"="#ffff99", "11"="#b15928") ###, "12"="#8dd3c7", "13"="#8e0152", "14"="#d9d9d9") 

mycol2 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
   "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
   "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")


p0 <- DimPlot(sc2, reduction="umap", cols=mycol2, group.by="predicted.id",
              label=T, repel=T, pt.size=0.2)+
    ## scale_color_manual(values=mycol2, drop=F)+
    theme_bw()+
    theme(plot.title=element_blank())

figfn <- "./3.2_annot.outs/Figure1.1_umap.celltype.png"
png(figfn, width=650, height=480, res=120)
p0
dev.off()


###
### heatmap

meta <- sc2@meta.data
meta <- meta%>%dplyr::select(NEW_BARCODE, seurat_clusters, predicted.id)

df1 <- meta%>%
   group_by(predicted.id, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))

p1 <- ggplot(df1, aes(x=seurat_clusters, y=predicted.id, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
                       low="#ffffc8", high="#7d0025", na.value=NA)+
   xlab("Cluster")+ylab("predicted.celltype")+
   theme_bw()+theme(axis.text.x=element_text(hjust=0.5, size=10))

###
png("./3.2_annot.outs/Figure1.2_heatmap.png", width=700, height=600, res=120)
print(p1)
dev.off()


###
### heatmap for cluster 17
sc3 <- sc2%>%FindClusters(resolution=0.1, verbose=T)

meta <- sc3@meta.data
meta <- meta%>%dplyr::select(NEW_BARCODE, seurat_clusters, predicted.id)

df1 <- meta%>%
   group_by(predicted.id, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))

p1 <- ggplot(df1, aes(x=seurat_clusters, y=predicted.id, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
                       low="#ffffc8", high="#7d0025", na.value=NA)+
   xlab("Cluster")+ylab("predicted.celltype")+
   theme_bw()+theme(axis.text.x=element_text(hjust=0.5, size=10))

###
png("./3.2_annot.outs/Figure1.3_heatmap_cluster17.png", width=800, height=600, res=120)
print(p1)
dev.off()



#############################
### second reference data ###
#############################
fn <- "/wsu/home/groups/bannonlab/sc_brains/Analysis/Ref2_nn2022/1_seurat.outs/1_seurat.rds"
ref <- read_rds(fn)
ref <- ref%>%
   NormalizeData()%>%
   FindVariableFeatures(nfeatures=2000)%>%
   ScaleData()%>%
   RunPCA(npcs=100)

opfn <- "/wsu/home/groups/bannonlab/sc_brains/Analysis/Ref2_nn2022/1_seurat.outs/2_seurat.rds"
write_rds(ref, opfn)

### transfer cell-type
anchors <- FindTransferAnchors(reference=ref,  query=sc, dims=1:50, reference.reduction="pca")
pred <- TransferData(anchorset=anchors,  refdata=ref$cellType2, dims=1:50)
pred$NEW_BARCODE <- rownames(pred)

opfn <- "./3.2_annot.outs/2_pred.rds"
write_rds(pred, opfn)


opfn <- "./3.2_annot.outs/2_anchors.rds"
write_rds(anchors, opfn)

#########################
### umap of cell-type ###
#########################

sc <- read_rds("./3_cluster.outs/1_seurat.cluster14.rds")
x <- sc@meta.data

pred <- read_rds("./3.2_annot.outs/2_pred.rds")
pred2 <- pred%>%dplyr::select(NEW_BARCODE, predicted.id)

x2 <- x%>%left_join(pred2, by="NEW_BARCODE")
rownames(x2) <- x2$NEW_BARCODE

sc2 <- AddMetaData(sc, metadata=x2)
Idents(sc2) <- sc2$predicted.id

## ##
## mycol1 <- c("0"="#8dd3c7", "1"="#ffffb3", "2"="#bebada", "3"="#fb8072",
##   "4"="#80b1d3", "5"="#fdb462", "6"="#b3de69", "7"="#fccde5", "8"="#d9d9d9",
##   "9"="#bc80bd", "10"="#ccebc5", "11"="#ffed6f")

## mycol2 <- c("0:ODC"="#a6cee3", "1:ODC"="#1f78b4", "2:ODC"="#b2df8a", "3:Astrocyte"="#33a02c",
##    "4:Microglia"="#fb9a99", "5:OPC"="#e31a1c", "6:Astrocyte"="#fdbf6f", "7:Dopamine Neurons"="#ff7f00",
##    "8:Endothelial"="#cab2d6", "9:GABA"="#6a3d9a")
##"10"="#ffff99", "11"="#b15928") ###, "12"="#8dd3c7", "13"="#8e0152", "14"="#d9d9d9") 

mycol2 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
   "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
   "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")


p0 <- DimPlot(sc2, reduction="umap", cols=mycol2, group.by="predicted.id",
              label=T, repel=T, pt.size=0.2)+
    ## scale_color_manual(values=mycol2, drop=F)+
    ggtitle("Automatic annotation")+
    theme_bw()+
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5))

figfn <- "./3.2_annot.outs/Figure2.1_umap.celltype.png"
png(figfn, width=440, height=480, res=120)
p0
dev.off()


###
### heatmap

meta <- sc2@meta.data
meta <- meta%>%dplyr::select(NEW_BARCODE, seurat_clusters, predicted.id)

df1 <- meta%>%
   group_by(predicted.id, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))

p1 <- ggplot(df1, aes(x=seurat_clusters, y=predicted.id, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
                       low="#ffffc8", high="#7d0025", na.value=NA)+
   xlab("Cluster")+ylab("predicted.celltype")+
   theme_bw()+theme(axis.text.x=element_text(hjust=0.5, size=10))

###
png("./3.2_annot.outs/Figure2.2_heatmap.png", width=700, height=600, res=120)
print(p1)
dev.off()


###
### heatmap for cluster 17
sc3 <- sc2%>%FindClusters(resolution=0.1, verbose=T)

meta <- sc3@meta.data
meta <- meta%>%dplyr::select(NEW_BARCODE, seurat_clusters, predicted.id)

df1 <- meta%>%
   group_by(predicted.id, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))

p1 <- ggplot(df1, aes(x=seurat_clusters, y=predicted.id, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
                       low="#ffffc8", high="#7d0025", na.value=NA)+
   xlab("Cluster")+ylab("predicted.celltype")+
   theme_bw()+theme(axis.text.x=element_text(hjust=0.5, size=10))

###
png("./3.2_annot.outs/Figure2.3_heatmap_cluster17.png", width=800, height=600, res=120)
print(p1)
dev.off()



######################################
### combine plots in supplementary ###
######################################

sc <- read_rds("./3_cluster.outs/1_seurat.cluster14.rds")
x <- sc@meta.data

pred <- read_rds("./3.2_annot.outs/2_pred.rds")
pred2 <- pred%>%dplyr::select(NEW_BARCODE, predicted.id)

x2 <- x%>%left_join(pred2, by="NEW_BARCODE")
rownames(x2) <- x2$NEW_BARCODE

sc2 <- AddMetaData(sc, metadata=x2)
Idents(sc2) <- sc2$predicted.id


###########################
###  umap cluster by ID ###
###########################

### final colors 
## col_MCls <- c("ODC"="#fb8072", "OPC"="#35978f", "Astrocyte"="#92cd2d", "Microglia"="#c38dc4", 
##    "DA"="#80b1d3", "Non_DA"="#fc9016",  "Pericyte"="#bebada",  "Endothelial"="#fccde5",
##     "T-cell"="#8dd3c7",  "Ependymal"="#828282")

###
### p1, umap colored by cluster 
mycol2 <- c("#fb8072", "#92cd2d", "#c38dc4", "#35978f", "#fc9016", "#b15928", "#bebada", "#80b1d3", "#fccde5",
   "#8dd3c7", "#ffff99", "#b2182b", "#828282", "#762a83")         

p1 <- DimPlot(sc, reduction="umap", cols=mycol2, group.by="seurat_clusters",
              label=T, repel=T, pt.size=0.2, label.size=3)+
    ##scale_color_manual(values=mycol2, drop=F)+
    ggtitle("Seurat clusters")+
    theme_bw()+
    theme(## legend.title=element_blank(),
          ## legend.key.size=grid::unit(0.5, "lines"),
          ## legend.text=element_text(size=10),
          legend.position="none",
          plot.title=element_text(hjust=0.5, size=12))


###
### p2, umap colored by cell-type from reference
mycol2 <- c("#92cd2d", "#80b1d3", "#fccde5", "#c38dc4", "#fc9016", "#fb8072", "#35978f")   

p2 <- DimPlot(sc2, reduction="umap", cols=mycol2, group.by="predicted.id",
              label=T, repel=T, pt.size=0.2, label.size=3)+
    ## scale_color_manual(values=mycol2, drop=F)+
    ggtitle("Automatic annotation")+
    theme_bw()+
    theme(## legend.title=element_blank(),
          ## legend.key.size=grid::unit(0.8, "lines"),
          ## legend.text=element_text(size=8),
          legend.position="none",
          plot.title=element_text(hjust=0.5, size=12))


###
### 3, heatmap show proportion

CL_value <- c("0"=1, "5"=2, "10"=3, "13"=4, "3"=5,
   "1"=6, "2"=7, "11"=8, "7"=9, "4"=10, "6"=11, "8"=12, "9"=13, "12"=14)

meta <- sc2@meta.data
meta <- meta%>%dplyr::select(NEW_BARCODE, seurat_clusters, predicted.id)

df1 <- meta%>%
   group_by(predicted.id, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))%>%ungroup()

df1 <- df1%>%
   mutate(CL_val=as.numeric(CL_value[as.character(seurat_clusters)]),
          Cluster2=fct_reorder(as.character(seurat_clusters), CL_val))          

p3 <- ggplot(df1, aes(x=Cluster2, y=predicted.id, fill=Perc))+
   geom_tile()+
   scale_fill_gradient(low="#ffffcc", high="#e31a1c",
     na.value=NA, limits=c(0,1), n.breaks=5,
     guide=guide_legend(keywidth=grid::unit(0.5,"lines"), keyheight=grid::unit(1,"lines") ))+
   xlab("Cluster id")+ ##ylab("Predicted.celltype")+
   ggtitle("Fraction of cells")+ 
   theme_bw()+
   theme(legend.title=element_blank(),
         axis.title.y=element_blank(),
         axis.text.x=element_text(hjust=0.5),
         plot.title=element_text(hjust=0.5, size=12))


###
### output

figfn <- paste(outdir, "Figure3_comb_ref2.png", sep="")
png(figfn, width=1100, height=400, res=120)
plot_grid(p1, p2, p3, ncol=3, labels="AUTO", label_fontface="plain", label_x=0.1,
          align="h", axis="tb", rel_widths=c(1,1,1.3))
dev.off()

