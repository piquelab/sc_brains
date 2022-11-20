##
library(Matrix)
library(tidyverse)
library(DESeq2)
library(annotables)
library(biobroom)
library(Seurat)
library(cowplot)
library(viridis)
##
rm(list=ls())

outdir <- "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

###
### Differential analysis based clusters, merging cluster 0, 5, 10 and 13 into one cluster 
### Last modified Aug-25-2022, by Julong wei

### Final model, model 2, ./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds
### model 3, ./4.0_Diff.cluster.outs/Filter2_cells30/1.2_DESeq.results_CPM0.5_filter.rds
### Filtering conditions, at least 30 cells per combination,
### autosome genes, CPM>0.5, at least 3 individuals at case and control respectively
### For each library, at least one case and control individual

########################
### pseudo-bulk data ###
########################

sc <- read_rds("./3_cluster.outs/1_seurat.cluster14.rds")

## counts data
## counts <- sc@assays$RNA@counts
## meta <- sc@meta.data
## cellSel <- meta%>%filter(seurat_clusters==0)%>%pull(NEW_BARCODE)
## count2 <- counts[,cellSel]
## rnz <- rowSums(count2)
## sum(rnz>20) ## try 0, 50, 100


anno <- data.frame(rn=rownames(counts), rnz=rowSums(counts))%>%filter(rnz>0)
    
## ##gene:64428, symbol:58149, ensgene:63697,
autosome <- as.character(1:22)
grch38_unq <- grch38%>%filter(symbol%in%anno$rn, chr%in%autosome)%>%
   distinct(symbol, chr, .keep_all=T)%>%
   dplyr::select(ensgene, symbol, chr, start, end, biotype) 
grch38_unq <- grch38_unq%>%group_by(symbol)%>%mutate(ny=n())%>%filter(ny==1)
 
annoSel <- anno%>%filter(rn%in%grch38_unq$symbol)  ### 28,327 genes
        
Y <- counts[annoSel$rn,]

## meta data
meta <- sc@meta.data%>%
   mutate(seurat_clusters=as.character(seurat_clusters),
          seurat_cluster2=ifelse(seurat_clusters%in%c("0", "5", "10", "13"), "0", seurat_clusters),
          seurat_cluster2=ifelse(seurat_cluster2%in%c("2", "11"), "2", seurat_cluster2))

meta <- meta%>%
       mutate(bti=paste(seurat_cluster2,  sampleID,  sep="_"))%>%
       dplyr::select(NEW_BARCODE, bti)

dd <- meta %>%group_by(bti)%>%summarise(ncell=n(),.groups="drop")

       
bti <- factor(meta$bti)       
X <- model.matrix(~0+bti)
YtX <- Y %*% X
YtX <- as.matrix(YtX)
colnames(YtX) <- gsub("^bti", "", colnames(YtX))
## opfn <- paste(outdir, "YtX.comb.rds", sep="")
## write_rds(YtX, opfn)

### Filter 30 cells per combination
ddx <- dd%>%filter(ncell>=30)
tmp <- ddx%>%mutate(MCls=gsub("_.*", "", bti), sampleID=gsub(".*_", "", bti))%>%as.data.frame()
tmp%>%group_by(MCls)%>%summarise(ncell=n(),.groups="drop")%>%arrange(as.numeric(MCls))
###
## YtX <- read_rds("./4.0_Diff.cluster.outs/Filter2_cells30/YtX.comb.rds")
## YtX_sel <- YtX[,ddx$bti]
opfn <- paste(outdir, "YtX_sel.comb.rds", sep="")
write_rds(YtX_sel, opfn)

### No filters
## opfn <- "./6_DEG.CelltypeNew_output/Filter2/YtX.comb.RData"
## save(YtX, file=opfn)
## ###


###
### covariance information
x <- openxlsx::read.xlsx("cell.counts_2022-07-06_final.xlsx")
x <- x[,c(1,3:17)]
names(x)[2] <- "sampleID"

### genotype PCs
pcs <- read.table("/wsu/home/groups/bannonlab/sc_brains/genotype_merge/plink/ref.eigenvec")
pcs <- pcs%>%mutate(sampleID=toupper(gsub(".*_", "", V2)))
pcs <- pcs%>%dplyr::select(sampleID, PC1=V3, PC2=V4, PC3=V5)

x2 <- x%>%left_join(pcs, by="sampleID")

write.csv(x2, "BrainCV_cell.counts.final.csv", row.names=F)


### number of cells in case and control 
x3 <- x2%>%filter(USE==1)%>%dplyr::select(sampleID, Category)
meta2 <- meta%>%left_join(x3, by="sampleID")
dd <- meta2%>%group_by(seurat_clusters, Category)%>%summarise(ny=n(),.groups="drop")
dd%>%pivot_wider(names_from=Category, values_from=ny, values_fill=0)%>%as.data.frame()

### data used for final data analysis
## fn <- "BrainCV_cell.counts.final.csv"
## cv <- read.csv(fn)%>%dplyr::filter(USE==1)

## x <- sc@meta.data
## x2 <- x%>%group_by(sampleID)%>%summarise(ncell=n(),.groups="drop")%>%as.data.frame()

## ncell <- x2$ncell
## names(ncell) <- as.character(x2$sampleID) 
## cv$n <- ncell[as.character(cv$sampleID)]

## opfn <- "BrainCV_cell_counts_final_USE_JW.xlsx"
## openxlsx::write.xlsx(cv, opfn, overwrite=T)




##############################
### runDESeq data analysis ###
##############################

fn <-  "./4.0_Diff.cluster.outs/YtX_sel.comb.rds"
YtX <- read_rds(fn)

### CPM 
depth <- colSums(YtX)
CPM <- sweep(YtX, 2, STATS=depth, FUN="/")*1e+06
## CPM_index <- CPM>0.5
 
## YtX2 <- YtX[rowSums(YtX)>20,]

bti2 <- colnames(YtX)
cvt0 <- str_split(bti2, "_", simplify=T)
cvt <- data.frame(bti=bti2, MCls=cvt0[,1], sampleID=cvt0[,2])
rownames(cvt) <- cvt$bti

##
x <- read.csv("BrainCV_cell.counts.final.csv")%>%dplyr::filter(USE==1)%>%
   dplyr::select(sampleID, Library, Age, Race, Sex, BW_lbs, BW_g, pH,
                 CauseOfDeath, BrainRegion, Category, PC1, PC2, PC3)
cvt <- cvt%>%left_join(x, by="sampleID")%>%
    mutate(Data_source=ifelse(grepl("AKB", sampleID), "Bannon", "Mash"))

## cvt_summ <- cvt%>%group_by(MCls, Category)%>%summarise(n=n(),.groups="drop")%>%
##   pivot_wider(names_from=Category, values_from=n)

## cvt_summ2 <- cvt%>%group_by(MCls, Library, Category)%>%summarise(n=n(),.groups="drop")%>%
##     mutate(comb=paste(MCls, Library, sep="_"))%>%
##     pivot_wider(id_cols=comb, names_from=Category, values_from=n, values_fill=0)%>%as.data.frame()


## cvt_summ2 <- cvt%>%filter(MCls=="6")%>%group_by(Library, Category)%>%summarise(n=n(),.groups="drop")%>%
##     pivot_wider(id_cols=Library, names_from=Category, values_from=n, values_fill=0)%>%as.data.frame()

## ###
## cvt_summ2 <- cvt_summ2%>%
##     mutate(MCls=gsub("_.*", "", comb), Library=gsub(".*_", "", comb))%>%
##     dplyr::select(MCls, Library, Control, Opioid)

## tmp <- cvt_summ2%>%filter(MCls==6) ## Control==0|Opioid==0)

## openxlsx::write.xlsx(cvt_summ2, "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/0_summary.xlsx", overwrite=T)


##
th_cpm <- 0.5
MCls <- as.character(sort(as.numeric(unique(cvt$MCls))))[1:9]
res <- lapply(MCls, function(ii){
   time0 <- Sys.time()
   ###
   cvt0 <- cvt%>%dplyr::filter(MCls==ii)  ## %>%mutate(Lib2=paste(Data_source, Library, sep="_"))
   ## lib_summ <- cvt0%>%group_by(Lib2)%>%summarize(ny=n(), .groups="drop")%>%as.data.frame()
   ## lib_summ <- lib_summ%>%
   ##     mutate(Lib2new=ifelse(ny==1, gsub("_.*", "", Lib2), gsub(".*_", "", Lib2)))%>%
   ##     dplyr::select(-ny)
   ##
   ## cvt0 <- cvt0%>%left_join(lib_summ, by="Lib2")

   lib_summ <- cvt0%>%group_by(Library, Category)%>%summarise(ny=n(), .groups="drop")%>%
       pivot_wider(id_cols=Library, names_from=Category, values_from=ny, values_fill=0)
   LibSel <- lib_summ%>%dplyr::filter(Control>=1, Opioid>=1)%>%dplyr::pull(Library)
   cvt0 <- cvt0%>%dplyr::filter(Library%in%LibSel)
       
   bti1 <- cvt0%>%filter(Category=="Control")%>%dplyr::pull(bti)
   bti2 <- cvt0%>%filter(Category=="Opioid")%>%dplyr::pull(bti)

   if ( length(bti1)>=3&length(bti2)>=3){
   
   geneSel <- rowSums(CPM[,bti1]>th_cpm)>=3&rowSums(CPM[,bti2]>th_cpm)>=3

   YtX0 <- YtX[geneSel, cvt0$bti]

   if ( length(unique(cvt0$Sex))>1){
       dds <- try(DESeqDataSetFromMatrix(YtX0, cvt0,
          ~Category+factor(Library)+as.numeric(Age)+as.numeric(pH)+as.numeric(PC1)+as.numeric(PC2)+as.numeric(PC3)+factor(Sex) ))
   }else{
       dds <- try(DESeqDataSetFromMatrix(YtX0, cvt0,
          ~Category+factor(Library)+as.numeric(Age)+as.numeric(pH)+as.numeric(PC1)+as.numeric(PC2)+as.numeric(PC3) ))
   }    
   dds <- try(DESeq(dds))
   if ( (class(dds)!="try-error") ){
      ###
      res <- results(dds, contrast=c("Category", "Opioid", "Control")) #independentFiltering=F
      res2 <- res%>%tidy()%>%mutate(MCls=ii, contrast="Opioid")
   }else{
      res2 <- NA
   }
   }else{
      res2 <- NA
   }
   ###
   time1 <- Sys.time()
   elapsed <- difftime(time1, time0, units="mins")
   cat(ii, elapsed, "Done\n")
   res2
})
####
###
res <- res[!is.na(res)]
res <- do.call(rbind,res)

##
opfn <- paste(outdir, "1.4_DESeq.results_CPM0.5_filter.rds", sep="")
write_rds(res, opfn)



################
### MA plots ###
################
 
MCls_label <- c("0"="ODC", "1"="astrocyte", "2"="microglia", "3"="OPC", "4"="GABA",
   "6"="pericyte", "7"="DaN", "8"="endothelial", "9"="T-cells")

###
### read data
fn <- paste(outdir, "1.2_DESeq.results_CPM0.5_filter.rds", sep="")
## fn <- paste(outdir, "1.0_DESeq.results.rds", sep="")
res <- read_rds(fn)%>%
   mutate(color=ifelse((p.adjusted<0.1)&(!is.na(p.adjusted)), T, F))


figfn <- paste(outdir, "Figure2.1_MA.png", sep="")
png(figfn, width=800, height=900, pointsize=12, res=150)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:9, 3, 3, byrow=T)
layout(x)

 
MCls <- as.character(sort(as.numeric(unique(res$MCls))))
for (ii in MCls){
   ##1
   oneMCl <- paste(ii, MCls_label[ii], sep=":")
   cat(oneMCl, "\n") 
   res2 <- res%>%
       dplyr::filter(MCls==as.character(ii))%>%
       dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main=oneMCl, cex.main=1.2, cex.axis=0.8, cex.lab=1))
}
dev.off()


###
### QQ plots
fn <- paste(outdir, "1.2_DESeq.results_CPM0.5_filter.rds", sep="")
res <- read_rds(fn)%>%drop_na(p.value)

###
MCls <- as.character(sort(as.numeric(unique(res$MCls))))
plotDF <- map_dfr(MCls, function(ii){
   ## 
   res2 <- res%>%filter(MCls==ii)%>%arrange(p.value)
   ngene <- nrow(res2)
   res2 <- res2%>%mutate(observed=-log10(p.value), expected=-log10(ppoints(ngene)))
   res2
})

plotDF <- plotDF%>%
    mutate(MCl_value=as.numeric(MCls),
    MCls=fct_reorder(MCls, MCl_value))
     
### qq plots
p  <- ggplot(plotDF, aes(x=expected, y=observed))+
   ggrastr::rasterise(geom_point(size=0.3, colour="grey30"),dpi=300)+
   facet_wrap(~MCls, ncol=3, scales="free",
     labeller=as_labeller(c("0"="0:ODC", "1"="1:astrocyte", "2"="2:microglia", "3"="3:OPC", "4"="4:GABA",
   "6"="6:pericyte", "7"="7:DaN", "8"="8:endothelial", "9"="9:T-cells")))+
   geom_abline(colour="red")+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(strip.text=element_text(size=12))
###

figfn <- paste(outdir, "Figure2.2_qq.png", sep="")
png(figfn, width=750, height=800, res=120)
print(p)
dev.off()



#####################
### Table of DEGs ###
#####################

outdir <- "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/"    
fn <- paste(outdir, "1.4_DESeq.results_CPM0.5_filter.rds", sep="")

## fn <- paste(outdir, "1.2_DESeq.results_CPM0.5_filter.rds", sep="")

res <- read_rds(fn)
 
### DEGs
res2 <- res%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.25)
opfn <- paste(outdir, "1.4_DESeq.xlsx", sep="")
openxlsx::write.xlsx(res2, opfn)


tmp <- res2%>%group_by(MCls)%>%summarise(ngene=n(), .groups="drop")%>%arrange(as.numeric(MCls))

length(unique(res2$gene))

## test genes
tmp <- res%>%group_by(MCls)%>%summarise(ngene=n(), .groups="drop")%>%arrange(as.numeric(MCls))
###
## x <- res%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.5)



###
###
res2 <- res%>%group_by(MCls)%>%arrange(desc(abs(estimate)), .by_group=T)
openxlsx::write.xlsx(res2, "./4.0_Diff.cluster.outs/1.1_DESeq.results.xlsx", overwrite=T)

###
res2 <- res%>%filter(abs(estimate)>0.25,p.adjusted<0.1)%>%
    group_by(MCls)%>%arrange(desc(abs(estimate)), .by_group=T)
openxlsx::write.xlsx(res2, "./4.0_Diff.cluster.outs/1.2_DESeq.LFC.0.25.xlsx", overwrite=T)


##################################################################
### scatter plots of z-score of sc-brain to previous bulk data ###
##################################################################

### function
feq <- function(x){
  r <- round(as.numeric(x$estimate),digits=3)
  p <- x$p.value
  if(p<0.001) symb <- "***"
  if(p>=0.001 & p<0.01) symb <- "**"
  if (p>=0.01 & p<0.05) symb <- "*"
  if(p>0.05) symb <- "NS"
  
  eq <- bquote(italic(R)==.(r)~","~.(symb))
  eq 
}

##
xFun <- function(dx,a=0.5){
min1 <- min(dx$x, na.rm=T)
max2 <- max(dx$x, na.rm=T)
R <- max2-min1
xpos <- min1+a*R
}
##
yFun <- function(dx, a=0.8){
min1 <- min(dx$y, na.rm=T)
max2 <- max(dx$y, na.rm=T)
R <- max2-min1
ypos <- min1+a*R
}



#####################################################
### set parameters that are required in the plots ###
#####################################################

### current color
col_MCls <- c("ODC"="#fb8072", "Astrocyte"="#92cd2d", "Microglia"="#c38dc4", "OPC"="#35978f",
   "Non_DA"="#fc9016",  "Pericyte"="#bebada", "DA"="#80b1d3", "Endothelial"="#fccde5",
    "T-cell"="#8dd3c7",  "Ependymal"="#828282")

### mapping cluster to cell-type
MCls_name <- c("0"="ODC", "1"="Astrocyte", "2"="Microglia", "3"="OPC",
     "4"="Non_DA", "6"="Pericyte", "7"="DA", "8"="Endothelial", "9"="T-cell", "12"="Ependymal")

### order of cell-type 
MCls_val <- c("ODC"=1, "OPC"=2, "Astrocyte"=3, "Microglia"=4, "DA"=5, "Non_DA"=6, 
    "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)



#################
### plot data ###
#################

## results from old bulk
fn <- "/wsu/home/groups/bannonlab/manal/DESeq2_WGCNA_DownstreamAnalysis/Final_Runs/DESEQ2_TEST_3-28-18/DESEQ2_IHW-5-28-18.csv"
old <- read.csv(fn)%>%drop_na(pvalue)%>%mutate(zscore_old=log2FoldChange/lfcSE)%>%    
   dplyr::select(gene=Associated.Gene.Name, beta_old=log2FoldChange, zscore_old, pval_old=pvalue) 

### results from sc
fn <- "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds"
res2 <- read_rds(fn)%>%drop_na(p.value)%>%
   mutate(zscore_sc=estimate/stderror)%>%
   dplyr::select(MCls, gene, beta_sc=estimate, zscore_sc, pval_sc=p.value)

plotDF <- res2%>%inner_join(old, by="gene")


###
### scatter plots of zscore between sc-brain and bulk

plotDF2 <- plotDF%>%dplyr::select(MCls, gene, x=zscore_sc, y=zscore_old)%>%
    mutate(celltype=MCls_name[as.character(MCls)],
           MCls_value=as.numeric(MCls_val[as.character(celltype)]),
           celltype2=fct_reorder(celltype, MCls_value))

###
### annotation data frame
anno_df1 <- plotDF2%>%group_by(celltype)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=map_dbl(data,~xFun(.x,a=0.5)),
          ypos=map_dbl(data,~yFun(.x,a=0.98)))%>%
   dplyr::select(-data,-corr)%>%unnest(celltype)

anno_df1 <- anno_df1%>%
    mutate(MCls_value=as.numeric(MCls_val[as.character(celltype)]),
           celltype2=fct_reorder(celltype, MCls_value))

###
### main figures
fig1 <- ggplot(plotDF2, aes(x=x, y=y))+ ## fill=factor(celltype)))+
   ##geom_point(size=0.3)+
   stat_density_2d(aes(fill=..level..), geom="polygon")+ 
   geom_text(data=anno_df1, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+
   scale_fill_viridis(option="C")+ 
   ##scale_fill_manual(values=col_MCls, guide="none")+
   facet_wrap(~celltype2, nrow=2, ncol=4, scales="free")+         
   scale_x_continuous("z-score from different cell-types in sc-brain", expand=expansion(mult=0.1))+
   scale_y_continuous("z-score from Bulk", expand=expansion(mult=0.15))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.position=c(0.9, 0.25),
         legend.key.size=grid::unit(0.8, "lines"),
         strip.text=element_text(size=12))
   
fig1 <- fig1+geom_smooth(method="lm",formula=y~x, size=0.5, se=F, color="red")
                           
figfn <- paste(outdir, "Figure4.4_score_SCvsBulk.png", sep="")
png(figfn, width=880, height=500, res=120)
print(fig1)
dev.off()

###
###

## fn <- "./4.0_Diff.cluster.outs/Filter2_cells30_M1_no/1.0_DESeq.results.rds"
## old <- read_rds(fn)%>%drop_na(p.value)%>%
##    mutate(zscore_old=estimate/stderror, MCls_gene=paste(MCls, gene, sep="_"))%>%    
##    dplyr::select(MCls_gene, zscore_old) 

## ### 
## outdir <- "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/"
## fn <- paste(outdir, "1.0_DESeq.results.rds", sep="")
## res2 <- read_rds(fn)%>%drop_na(p.value)%>%
##    mutate(zscore_sc=estimate/stderror, MCls_gene=paste(MCls, gene, sep="_"))%>%
##    dplyr::select(MCls, gene,  MCls_gene, zscore_sc)

## plotDF <- res2%>%inner_join(old, by="MCls_gene")


## ###
## ### scatter plots of zscore between sc-brain and bulk

## plotDF2 <- plotDF%>%dplyr::select(MCls, gene, MCls_gene, x=zscore_old, y=zscore_sc)%>%
##     mutate(MCls_value=as.numeric(MCls),
##            MCls_label=paste("Cluster", MCls),
##            MCls_label=fct_reorder(MCls_label, MCls_value))
 
## anno_df1 <- plotDF2%>%group_by(MCls)%>%
##    nest()%>%
##    mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
##           eq=map(corr,feq),
##           r2=map_dbl(corr,~(.x)$estimate),
##           xpos=map_dbl(data,~xFun(.x,a=0.5)),
##           ypos=map_dbl(data,~yFun(.x,a=0.98)))%>%
##    dplyr::select(-data,-corr)

## anno_df1 <- anno_df1%>%
##     mutate(MCls_value=as.numeric(MCls),
##            MCls_label=paste("Cluster", MCls),
##            MCls_label=fct_reorder(MCls_label, MCls_value))

## fig1 <- ggplot(plotDF2, aes(x=x, y=y))+
##    geom_point(size=0.3, color="grey50")+ 
##    geom_text(data=anno_df1, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+ 
##    facet_wrap(~MCls_label, ncol=3, scales="free")+         
##    scale_x_continuous("zscore of model w/o library", expand=expansion(mult=0.1))+
##    scale_y_continuous("zscore from w/. library", expand=expansion(mult=0.15))+
##    theme_bw()+
##    theme(strip.text=element_text(size=12),
##          axis.title=element_text(size=12))
## fig1 <- fig1+geom_smooth(method="lm",formula=y~x, size=0.5, se=F, color="red")
                           
## figfn <- paste(outdir, "Figure2.1_score_Model_2vs1.png", sep="")
## png(filename=figfn, width=750, height=900, res=120)  
## print(fig1)
## dev.off()


###
###
fn <-  "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds"
res <- read_rds(fn)

summ <- res%>%group_by(MCls)%>%
   summarise(ntest=sum(!is.na(estimate)),
             ntest_q=sum(!is.na(p.adjusted)), n_na=sum(is.na(p.adjusted)), .groups="drop")%>%
   ungroup()

opfn <- "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/2_summary_qc.xlsx"
openxlsx::write.xlsx(summ, file=opfn, overwrite=T)
### scatter plots of log2FoldChange

## plotDF2 <- plotDF%>%dplyr::select(MCls, gene, x=beta_sc, y=beta_old)%>%
##     mutate(MCls_value=as.numeric(MCls),
##            MCls_label=paste("Cluster", MCls),
##            MCls_label=fct_reorder(MCls_label, MCls_value))
 
## anno_df1 <- plotDF2%>%group_by(MCls)%>%
##    nest()%>%
##    mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
##           eq=map(corr,feq),
##           r2=map_dbl(corr,~(.x)$estimate),
##           xpos=map_dbl(data,~xFun(.x,a=0.5)),
##           ypos=map_dbl(data,~yFun(.x,a=0.98)))%>%
##    dplyr::select(-data,-corr)

## anno_df1 <- anno_df1%>%
##     mutate(MCls_value=as.numeric(MCls),
##            MCls_label=paste("Cluster", MCls),
##            MCls_label=fct_reorder(MCls_label, MCls_value))

## fig1 <- ggplot(plotDF2, aes(x=x, y=y))+
##    geom_point(size=0.3, color="grey50")+ 
##    geom_text(data=anno_df1, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+ 
##    facet_wrap(~MCls_label, ncol=3, scales="free")+         
##    scale_x_continuous("log2FoldChange from sc-brain", expand=expansion(mult=0.1))+
##    scale_y_continuous("log2FoldChange from bulk", expand=expansion(mult=0.15))+
##    theme_bw()+
##    theme(strip.text=element_text(size=12),
##          axis.title=element_text(size=12))
## fig1 <- fig1+geom_smooth(method="lm", formula=y~x, size=0.5, se=F, color="red")
                           
## figfn <- paste(outdir, "Figure2.0_log2FoldChange_bulk.png", sep="")
## png(filename=figfn, width=750, height=900, res=120)  
## print(fig1)
## dev.off()


###
### scatter plots of pvalue between sc-brain and previous bulk

## plotDF2 <- plotDF%>%dplyr::select(MCls, gene, x=pval_sc, y=pval_old)%>%
##     mutate(MCls_value=as.numeric(MCls),
##            MCls_label=paste("Cluster", MCls),
##            MCls_label=fct_reorder(MCls_label, MCls_value))

## fig1 <- ggplot(plotDF2, aes(x=-log10(x), y=-log10(y)))+
##    geom_point(size=0.3, color="grey50")+ 
##    geom_abline(slope=1, intercept=0, color="red")+
##    facet_wrap(~MCls_label, ncol=3, scales="free")+         
##    scale_x_continuous(bquote(~Log[10]~" pvalue from sc-brain"), expand=expansion(mult=0.1))+
##    scale_y_continuous(bquote(~Log[10]~" pvalue from bulk"), expand=expansion(mult=0.15))+
##    theme_bw()+
##    theme(strip.text=element_text(size=12),
##          axis.title=element_text(size=12))

                           
## figfn <- paste(outdir, "Figure2.2_pval_bulk.png", sep="")
## png(filename=figfn, width=750, height=900, res=120)  
## print(fig1)
## dev.off()


################################
#### forest plots of 10 DEGs ###
################################
 
## fn <- paste(outdir, "1.4_DESeq.results_CPM0.5_filter.rds", sep="")
## res <- read_rds(fn)%>%drop_na(p.value)
 
## DEGlist <- c("PLA1A", "MAP3K6", "YBX3", "FOSL2", "MGP", "MEFV", "CISH", "PDLIM1", "GPR4", "TRIP10")
## res2 <- res%>%filter(gene%in%DEGlist)%>%
##    dplyr::select(gene, beta=estimate, SE=stderror, MCls, p.adjusted)


## mycol2 <- c("0"="#a6cee3", "1"="#1f78b4", "2"="#b2df8a", "3"="#33a02c",
##    "4"="#fb9a99", "5"="#e31a1c", "6"="#fdbf6f", "7"="#ff7f00", "8"="#cab2d6",
##     "9"="#6a3d9a","10"="#ffff99", "11"="#b15928", "not_DEG"="grey40") #, "12"="#8dd3c7", "13"="#8e0152", "14"="#d9d9d9")
 
## plotDF <- res2%>%
##    mutate(beta_upper=beta+1.96*SE, beta_lower=beta-1.96*SE,
##           gr2=ifelse(p.adjusted<0.1&abs(beta)>0.25, as.character(MCls), "not_DEG"))

## p <- ggplot(plotDF, aes(x=beta, y=factor(MCls), color=factor(gr2)))+
##    geom_errorbarh(aes(xmax=beta_upper, xmin=beta_lower), size=0.5, height=0.2)+
##    geom_point(shape=19, size=1.5)+
##    scale_colour_manual(values=mycol2)+
##    xlab("Log2FoldChange")+
##    scale_y_discrete("",
##       labels=c("0"="0:ODC", "1"="1:astrocyte", "2"="2:microglia", "3"="3:OPC", "4"="4:GABA",
##                "6"="6:pericyte", "7"="7:DaN", "8"="8:endothelial", "9"="9:T-cells"))+
##    geom_vline(xintercept=0, size=0.25, linetype="dashed")+
##    facet_wrap(~factor(gene), ncol=5, scales="free")+
##    theme_bw()+
##    theme(legend.position="none")
                           
## figfn <- paste(outdir, "Figure4.5_forestplot.png", sep="")
## png(filename=figfn, width=1200, height=600, res=120)  
## print(p)
## dev.off()
   


##############################
#### Plots for publication ###
##############################


## outdir <- "./4.0_Diff.cluster.outs/plots_pub/"
## if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

## fn <-  "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds"
## res <- read_rds(fn)
 

## col_MCls <- c("ODC"="#fb8072", "Astrocyte"="#b3de69", "Microglia"="#bebada", "OPC"="#8dd3c7",
##    "Non_DA"="#fdb462",  "Pericyte"="#ffffb3", "DA"="#80b1d3", "Endothelial"="#fccde5",
##     "T-cell"="#35978f", "Ependymal"="#828282")
## MCls_name <- c("0"="ODC", "1"="Astrocyte", "2"="Microglia", "3"="OPC",
##      "4"="Non_DA", "6"="Pericyte", "7"="DA", "8"="Endothelial", "9"="T-cell", "12"="Ependymal")
## MCls_val <- c("ODC"=0, "Astrocyte"=1, "Microglia"=2, "OPC"=3, "Non_DA"=4, "Pericyte"=6)

## ### DEGs
## res2 <- res%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.25)%>%
##    mutate(direction=ifelse(estimate>0, 2, 1), celltype=MCls_name[as.character(MCls)])


## sigs <- res2%>%group_by(celltype, direction)%>%
##        summarise(ny=n(), .groups="drop")%>%
##        mutate(ny2=ifelse(direction==1, -ny, ny))%>%as.data.frame()
## sigs <- sigs%>%mutate(MCls_value=as.numeric(MCls_val[celltype]),
##                       celltype2=fct_reorder(celltype, MCls_value))

## breaks_value <- pretty(c(-1100, 1100), 5)

## p <- ggplot(sigs, aes(x=celltype2, y=ny2))+
##    geom_bar(aes(fill=factor(celltype), alpha=factor(direction)), stat="identity")+
##    scale_fill_manual(values=col_MCls)+
##    scale_alpha_manual(values=c("1"=0.5, "2"=1))+
##    geom_hline(yintercept=0, color="grey60")+
##    geom_text(aes(x=celltype, y=ny2, label=abs(ny2),
##        vjust=ifelse(direction==2, -0.2, 1.2)), size=2.8)+
##    scale_y_continuous("", breaks=breaks_value, limits=c(-1150,1150), labels=abs(breaks_value))+
##    theme_bw()+
##    theme(legend.position="none",
##          axis.title=element_blank(),
##          axis.text=element_text(size=10))

## ###
## figfn <- paste(outdir, "Figure1.1_barplot.png", sep="")
## png(filename=figfn, width=580, height=400, res=120)
## print(p)
## dev.off()







################################################################
### mapping reads to X/Y chromsome and check Sex information ###
################################################################

## sc <- read_rds("./2_seurat.outs/2_seurat.clean.rds")

## ## counts data 
## Y <- sc@assays$RNA@counts


## ## filtering genes
## anno <- data.frame(rn=rownames(Y), rnz=rowSums(Y))
    
## ## 1,202 gene
## grch38_unq <- grch38%>%filter(symbol%in%anno$rn, chr%in%c("X", "Y"))%>%
##    distinct(symbol, chr, .keep_all=T)%>%
##    dplyr::select(ensgene, symbol, chr, start, end, biotype)

## grch38_unq <- grch38_unq%>%group_by(symbol)%>%mutate(ny=n())%>%filter(ny==1)
## grch2 <- grch38_unq%>%dplyr::select(symbol, chr)

## annoSel <- anno%>%inner_join(grch2, by=c("rn"="symbol"))  ### 1,202 genes

## Ysel <- Y[annoSel$rn,]
## rownames(Ysel) <- paste(annoSel$chr, annoSel$rn, sep="_")


## ## meta data
## meta <- sc@meta.data
## bti <- factor(meta$sampleID)
## ### pseudo-bulk counts
## X <- model.matrix(~0+bti)
## YtX <- Ysel %*% X
## YtX <- as.matrix(YtX)
## colnames(YtX) <- gsub("^bti", "", colnames(YtX))

## ## opfn <- "./4.0_Diff.cluster.outs/YtX_chrXY.rds"
## ## write_rds(YtX, opfn)


## ## ###
## chr <- as.factor(annoSel$chr)
## xchr <- model.matrix(~0+chr)

## reads <- t(YtX)%*%xchr

## Df <- reads%>%as.data.frame()%>%rownames_to_column(var="sampleID")


## ###
## ## dd <- meta %>%group_by(sampleID)%>%summarise(ncell=n(),.groups="drop")

## ## 
## x <- read.csv("BrainCV_cell.counts.csv")%>%filter(USE==1)%>%
##    dplyr::select(sampleID=Sample_ID2, Library,  Sex, Category)
 
## Df2 <- Df%>%left_join(x, by="sampleID")

## ##
## p1 <- ggplot(Df2)+
##    geom_point(aes(x=factor(sampleID), y=chrY, color=factor(Sex)))+
##    ylab("Reads from Y chromosome")+
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          axis.title.x=element_blank(),
##          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=6))
## ##
## figfn <- "./4.0_Diff.cluster.outs/Figure0.1_mapping_chrY.png"
## png(figfn, width=950, height=400, res=100)
## print(p1)
## dev.off()

## ##
## ##
## p2 <- ggplot(Df2)+
##    geom_point(aes(x=factor(sampleID), y=chrX, color=factor(Sex)))+
##    ylab("Reads from X chromosome")+
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          axis.title.x=element_blank(),
##          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=6))
## ##
## figfn <- "./4.0_Diff.cluster.outs/Figure0.2_mapping_chrX.png"
## png(figfn, width=950, height=400, res=100)
## print(p2)
## dev.off()


## ### scale 
## ### Need scale
## YtX_all <- Y%*%X
## depth <- colSums(YtX_all)
## depth <- data.frame(sampleID=gsub("^bti", "", names(depth)), depth=depth)

## Df3 <- Df2%>%left_join(depth, by="sampleID")%>%
##     mutate(chrX_2=(chrX/depth)*1e+04, chrY_2=(chrY/depth)*1e+04)

## p1 <- ggplot(Df3)+
##    geom_point(aes(x=factor(sampleID), y=chrY_2, color=factor(Sex)))+
##    ylab("Reads from Y chromosome (CP10k)")+
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          axis.title.x=element_blank(),
##          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=6))
## ##
## figfn <- "./4.0_Diff.cluster.outs/Figure0.1.1_scale_chrY.png"
## png(figfn, width=950, height=400, res=100)
## print(p1)
## dev.off()

## ##
## ##
## p2 <- ggplot(Df3)+
##    geom_point(aes(x=factor(sampleID), y=chrX_2, color=factor(Sex)))+
##    ylab("Reads from X chromosome (CP10k)")+
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          axis.title.x=element_blank(),
##          axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=6))
## ##
## figfn <- "./4.0_Diff.cluster.outs/Figure0.2.1_scale_chrX.png"
## png(figfn, width=950, height=400, res=100)
## print(p2)
## dev.off()



   
   
