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
library(DESeq2)
library(annotables)
##
library(ggrastr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(openxlsx)

## theme_set(theme_grey())

rm(list=ls())

#######################################################
### 1. finer cluster analysis for ODC and Microglia ###
#######################################################
 
outdir <- "./3_re-analysis.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

###
###
sc <- read_rds("../3_cluster.outs/1.2_seurat.annot.rds")


###
### re-analyse ODC
sc2 <- subset(sc, subset=celltype=="ODC")
x <- sc2@assays$RNA@counts
meta <- sc2@meta.data

###
sc3 <- CreateSeuratObject(x, project="sc-odc", meta.data=meta)

###
### normalized data, dimensional reduction and cluster analysis
sc3 <- sc3%>%
   NormalizeData()%>%
   FindVariableFeatures(nfeatures=2000)%>%
   ScaleData()%>%
   RunPCA(npcs=100)

sc3 <- RunHarmony(sc3, group.by.vars="EXP")
## sc <- sc%>%
##    RunUMAP(reduction="harmony", dims=1:50)%>%
##    FindNeighbors(dims=1:50, verbose=T)
sc3 <- sc3%>%
   RunUMAP(reduction="harmony", dims=1:50)%>%
   FindNeighbors(reduction="harmony", dims=1:50, verbose=T)

###
sc_final <- sc3%>%FindClusters(resolution=0.05, verbose=T)
 
opfn <- paste(outdir, "1.1_ODC_seurat.cluster.rds", sep="")
write_rds(sc_final, opfn)

###
### visulization data
p <- DimPlot(sc_final, reduction="umap", label=T)+theme_bw()
figfn <- paste(outdir, "Figure1.1_cl3_ODC_finer_umap.png", sep="")
png(figfn, width=520, height=480, res=120)
p
dev.off()


###
### subset Microglia and re-analysis

sc2 <- subset(sc, subset=celltype=="Microglia")
x <- sc2@assays$RNA@counts
meta <- sc2@meta.data
 
###
sc3 <- CreateSeuratObject(x, project="sc-microglia", meta.data=meta)

###
### rountine analysis
sc3 <- sc3%>%
   NormalizeData()%>%
   FindVariableFeatures(nfeatures=2000)%>%
   ScaleData()%>%
   RunPCA(npcs=100)

sc3 <- RunHarmony(sc3, group.by.vars="EXP")

sc3 <- sc3%>%
   RunUMAP(reduction="harmony", dims=1:50)%>%
   FindNeighbors(reduction="harmony", dims=1:50, verbose=T)

###
sc_final <- sc3%>%FindClusters(resolution=0.05, verbose=T)

### output
opfn <- paste(outdir, "1.2_Microglia_seurat.cluster.rds", sep="")
write_rds(sc_final, opfn)

### visulization data
p <- DimPlot(sc_final, reduction="umap", label=T)+theme_bw()
figfn <- paste(outdir, "Figure1.2_cl6_Microglia_finer_umap.png", sep="")
png(figfn, width=520, height=480, res=120)
p
dev.off()




#######################################
### Find cluster differential genes ###
#######################################

### ODC subclusters
fn <- "./3_re-analysis.outs/1.1_ODC_seurat.cluster.rds"
sc <- read_rds(fn)

resDiff <- FindAllMarkers(sc, only.pos=T, min.pct=0.02, logfc.threshold=0.1)

opfn <- "./3_re-analysis.outs/3.1_ODC_marker_genes.rds"
write_rds(resDiff, opfn)

###
### excel output
fn <- "./3_re-analysis.outs/3.1_ODC_marker_genes.rds"
res2 <- read_rds(fn)%>%filter(cluster%in%c(0,1))%>%
    group_by(cluster)%>%
    slice_max(order_by=avg_log2FC, n=100)

opfn <- paste(outdir, "3.1_ODC_marker_top100.xlsx", sep="")
write.xlsx(res2, opfn)


###
### microglia subclusters
fn <- "./3_re-analysis.outs/1.2_Microglia_seurat.cluster.rds"
sc <- read_rds(fn)

resDiff <- FindMarkers(sc, ident.1=0, ident.2=2, min.pct=0.02, logfc.threshold=0.1)

opfn <- "./3_re-analysis.outs/3.2_Microglia_marker_genes.rds"
write_rds(resDiff, opfn)

### excel output
fn <- "./3_re-analysis.outs/3.2_Microglia_marker_genes.rds"
res <- read_rds(fn)
res2 <- res%>%
    mutate(cluster=ifelse(sign(avg_log2FC)>0, "0", "2"), avg_log2FC=abs(avg_log2FC), gene=rownames(res))

res2 <- res2%>%group_by(cluster)%>%
    slice_max(order_by=avg_log2FC, n=100)%>%ungroup()%>%
    arrange(cluster)

opfn <- paste(outdir, "3.2_Microglia_marker_top100.xlsx", sep="")
write.xlsx(res2, opfn, overwrite=T)





#####################################
### 2. Generate pseudo-bulk data ####
#####################################

rm(list=ls())

outdir <- "./3_re-analysis.outs/Diff_analysis.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

sc <- read_rds("./3_re-analysis.outs/1.1_ODC_seurat.cluster.rds")

counts <- sc@assays$RNA@counts
meta <- sc@meta.data

anno <- data.frame(rn=rownames(counts), rnz=rowSums(counts))%>%filter(rnz>0)
    
## ##gene:64428, symbol:58149, ensgene:63697,
autosome <- as.character(1:22)
grch38_unq <- grch38%>%dplyr::filter(symbol%in%anno$rn, chr%in%autosome)%>%
   distinct(symbol, chr, .keep_all=T)%>%
   dplyr::select(ensgene, symbol, chr, start, end, biotype) 
grch38_unq <- grch38_unq%>%group_by(symbol)%>%mutate(ny=n())%>%filter(ny==1)
 
annoSel <- anno%>%filter(rn%in%grch38_unq$symbol)  ### 29,817
        
Y <- counts[annoSel$rn,]

## meta data
cluster_sel <- c("0", "1")
meta <- sc@meta.data%>%
   mutate(seurat_cluster2=as.character(seurat_clusters))%>%filter(seurat_cluster2%in%cluster_sel)

meta <- meta%>%
       mutate(bti=paste(seurat_cluster2,  sampleID,  sep="_"))%>%
       dplyr::select(NEW_BARCODE, bti)

dd <- meta %>%group_by(bti)%>%summarise(ncell=n(),.groups="drop")

Y <- Y[, meta$NEW_BARCODE]       
bti <- factor(meta$bti)       
X <- model.matrix(~0+bti)
YtX <- Y %*% X
YtX <- as.matrix(YtX)
colnames(YtX) <- gsub("^bti", "", colnames(YtX))
## opfn <- paste(outdir, "YtX.comb.rds", sep="")
## write_rds(YtX, opfn)


###
### Filter 30 cells per combination
ddx <- dd%>%filter(ncell>=30)
YtX_sel <- YtX[,ddx$bti]
opfn <- "./3_re-analysis.outs/2.1_ODC_YtX_sel.comb.rds"
write_rds(YtX_sel, opfn)



#################
### Microglia ### 
#################

sc <- read_rds("./3_re-analysis.outs/1.2_Microglia_seurat.cluster.rds")

counts <- sc@assays$RNA@counts
meta <- sc@meta.data

anno <- data.frame(rn=rownames(counts), rnz=rowSums(counts))%>%filter(rnz>0)
    
## ##gene:64428, symbol:58149, ensgene:63697,
autosome <- as.character(1:22)
grch38_unq <- grch38%>%dplyr::filter(symbol%in%anno$rn, chr%in%autosome)%>%
   distinct(symbol, chr, .keep_all=T)%>%
   dplyr::select(ensgene, symbol, chr, start, end, biotype) 
grch38_unq <- grch38_unq%>%group_by(symbol)%>%mutate(ny=n())%>%filter(ny==1)
 
annoSel <- anno%>%filter(rn%in%grch38_unq$symbol)  ### 28,623
        
Y <- counts[annoSel$rn,]

###
### meta data
meta <- sc@meta.data%>%
   mutate(seurat_cluster2=as.character(seurat_clusters),
          seurat_cluster2=ifelse(seurat_cluster2%in%c("0", "5"), "0", seurat_cluster2),
          seurat_cluster2=ifelse(seurat_cluster2%in%c("3", "4"), "3", seurat_cluster2))

meta <- meta%>%
       mutate(bti=paste(seurat_cluster2,  sampleID,  sep="_"))%>%
       dplyr::select(NEW_BARCODE, bti)

dd <- meta %>%group_by(bti)%>%summarise(ncell=n(),.groups="drop")

Y <- Y[, meta$NEW_BARCODE]       
bti <- factor(meta$bti)       
X <- model.matrix(~0+bti)
YtX <- Y %*% X
YtX <- as.matrix(YtX)
colnames(YtX) <- gsub("^bti", "", colnames(YtX))
## opfn <- paste(outdir, "YtX.comb.rds", sep="")
## write_rds(YtX, opfn)


###
### Filter 30 cells per combination
ddx <- dd%>%filter(ncell>=30)
YtX_sel <- YtX[,ddx$bti]
opfn <- "./3_re-analysis.outs/2.2_Microglia_YtX_sel.comb.rds"
write_rds(YtX_sel, opfn) 


