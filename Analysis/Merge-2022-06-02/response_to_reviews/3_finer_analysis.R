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

outdir <- "./3_re-analysis.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



#######################################################
### 1. finer cluster analysis for ODC and Microglia ###
#######################################################
 
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



####################################
## overlap marker gene expression ##
####################################
fn <- "./3_re-analysis.outs/3.1_ODC_marker_top100.xlsx"
x <- read.xlsx(fn)

###
fn <- "./3_re-analysis.outs/marques_ODC_aaf6463_table_s1.xlsx"
ref <- read.xlsx(fn)
ref2 <- ref[2:51,1:24]
colnames(ref2) <- ref[1,1:24]

###
MCls <- gsub(" |/|-", "_", colnames(ref2))
n_MCl <- length(MCls)


DF <- NULL
for (ii in c("0", "1")){
   ##
   x0 <- x%>%filter(cluster==ii)%>%dplyr::select(-p_val, -p_val_adj, -cluster)
   gene0 <- x0$gene 
   ##
   for (k in 1:n_MCl){
       ##
       gene1 <- toupper(ref2[,k])
       olap <- intersect(gene0, gene1)
       ##
       if ( length(olap)>0){
           ##
           df <- data.frame(cluster=ii, celltype=MCls[k], olap_gene=olap)
           df <- df%>%left_join(x0, by=c("olap_gene"="gene"))           
           DF <- rbind(DF, df)
       }    
   } ##
}

opfn <- paste(outdir, "3.1_ODC_olap.xlsx", sep="")
write.xlsx(DF, opfn, overwrite=T)
###

summ <- DF%>%group_by(cluster, celltype)%>%summarize(ngene=n(), .groups="drop")%>%ungroup()
summ2 <- summ%>%group_by(cluster)%>%arrange(desc(ngene), .by_group=T)
opfn <- paste(outdir, "3.1_ODC_olap_ngene.xlsx", sep="")
write.xlsx(summ2, opfn)


###
### Heatmap
fn <- "./3_re-analysis.outs/3.1_ODC_olap.xlsx"
DF <- read.xlsx(fn)

DF2 <- DF%>%group_by(cluster)%>%dplyr::distinct(olap_gene, .keep_all=T)%>%ungroup()
DF2 <- DF2%>%group_by(cluster)%>%arrange(desc(avg_log2FC), .by_group=T)%>%ungroup()

geneSel <- unique(DF2$olap_gene)

fn <- "./3_re-analysis.outs/1.1_ODC_seurat.cluster.rds"
sc <- read_rds(fn)
sc0 <- subset(sc, subset=seurat_clusters%in%c(0, 1))


p0 <- DoHeatmap(sc0, features=geneSel, slot="data", angle=0, size=1.5)+
   guides(color="none")+ 
   theme(axis.text=element_text(size=4),
         legend.key.size=grid::unit(0.2, "cm"),
         legend.text=element_text(size=5),
         legend.title=element_text(size=5))
 
figfn <- paste(outdir, "Figure3_ODC_heatmap.png", sep="")
ggsave(filename=figfn, plot=p0, device="png", width=8.5, height=4.8, units="cm")



####
###
fn <- "./3_re-analysis.outs/1.1_ODC_seurat.cluster.rds"
sc <- read_rds(fn)

x <- sc@meta.data
umap <- Embeddings(sc, reduction="umap")

## covariates
fn <- "../BrainCV_cell_counts_final_correct_2022-12-06.xlsx"
cvt <- read.xlsx(fn)
cvt2 <- cvt%>%dplyr::select(sampleID, Category)

plotDF <- data.frame(umap[,1:2], cluster=as.character(x$seurat_clusters), sampleID=x$sampleID)
plotDF <- plotDF%>%left_join(cvt2, by="sampleID")


##
### umap plots
lab_facet <- as_labeller(c("0"="cluster 0", "1"="cluster 1"))
plotDF2 <- plotDF%>%filter(cluster%in%c("0", "1"))
###
p <- ggplot(plotDF2, aes(x=UMAP_1, y=UMAP_2))+
   geom_point(aes(color=factor(Category)), size=0.1)+
   scale_colour_manual(values=c("Control"="#0571b0", "Opioid"="#ca0020"),
       guide=guide_legend(override.aes=list(size=3)))+
   facet_wrap(~factor(cluster), labeller=lab_facet, scales="fixed", ncol=2)+ 
   theme_bw()+
   theme(axis.text=element_text(size=12),
         axis.title=element_text(size=12),
         legend.title=element_blank(),
         legend.key.size=unit(0.5, "cm"),
         strip.text=element_text(size=14))

figfn <- paste(outdir, "Figure3.1_subcluster.png", sep="")
png(figfn, width=680, height=380, res=120)
print(p)
dev.off()

summ <- plotDF2%>%group_by(cluster, Category)%>%summarize(ncell=n(),.groups="drop")


################################ 
### OPALIN, S100B and RBFOX1 ###
################################

fn <- "./3_re-analysis.outs/Diff_analysis.outs/0.1_ODC_YtX_sel.comb.rds"
x <- read_rds(fn)
depth <- colSums(x)
x2 <- sweep(x, 2, depth, "/")*1e+6


###
fn <- "../BrainCV_cell_counts_final_correct_2022-12-06.xlsx"
cvt <- read.xlsx(fn)
cvt2 <- cvt%>%dplyr::select(sampleID, Category)



tmp <- str_split(colnames(x2), "_", simplify=T)
plotDF <- data.frame(rn=colnames(x2), MCls=tmp[,1], sampleID=tmp[,2]) 
plotDF <- plotDF%>%left_join(cvt2, by="sampleID")



geneList <- c("OPALIN", "S100B", "RBFOX1")
ngene <- length(geneList)
figs_ls <- lapply(1:ngene, function(i){
    ##
    gene <- geneList[i]
    cat(gene, "\n")
    plotDF2 <- plotDF%>%mutate(y=log2(x2[gene,]+1))
    ###
    p <- ggplot(plotDF2, aes(x=MCls, y=y, color=Category))+
       geom_violin(width=0.5)+
       geom_jitter(width=0.2, size=0.1)+
       ylab(bquote(~"Normalized gene expression"~"("~log[2]~"CPM"~")"))+
       scale_color_manual(values=c("Control"="#0571b0", "Opioid"="#ca0020"))+ 
       ##scale_fill_manual(values=c("Control"="#0571b0", "Opioid"="#ca0020"))+ 
       ggtitle(gene)+
       theme_bw()+
       theme(plot.title=element_text(hjust=0.5, size=12),
             axis.text.x=element_text(size=10),
             axis.text.y=element_text(size=10),
             axis.title.x=element_blank(),
             axis.title.y=element_text(size=10),
             legend.title=element_blank())
    ##
    if ( i==3){
       p <- p+theme(legend.key.size=unit(0.2, "cm"))
    }else{
       p <- p+theme(legend.position="none")
    }
    p
})


### output figures
figfn <- paste(outdir, "Figure4_marker_ODC_indi_violin.png", sep="")
png(figfn, width=780, height=380, res=120)
plot_grid(plotlist=figs_ls, ncol=3, rel_widths=c(1, 1, 1.4))
dev.off()


###
###
fn <- "./3_re-analysis.outs/Diff_analysis.outs/1.1_ODC_DESeq_results_default.rds"
res <- read_rds(fn)
res2 <- res%>%filter(abs(estimate)>0.25, p.adjusted<0.1)

