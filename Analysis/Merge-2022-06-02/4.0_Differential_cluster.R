##
library(Matrix)
library(tidyverse)
library(DESeq2)
library(annotables)
library(biobroom)
library(Seurat)
library(cowplot)
library(viridis)
library(openxlsx)
##
rm(list=ls())


###
### Differential analysis based clusters, merging cluster 0, 5, 10 and 13 into one cluster 

###
### analysis-1st
### Last modified Aug-25-2022, by Julong wei
### Final model, model 2, ./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds
### model 3, ./4.0_Diff.cluster.outs/Filter2_cells30/1.2_DESeq.results_CPM0.5_filter.rds
### Filtering conditions, at least 30 cells per combination,
### autosome genes, CPM>0.5, at least 3 individuals at case and control respectively
### For each library, at least one case and control individual

## outdir <- "./4.0_Diff.cluster.outs/Old_results/Filter2_cells30_M2_batch/"
## fn <- paste(outdir, "1.4_DESeq.results_CPM0.5_filter.rds", sep="")
## res <- read_rds(fn)

###
### analysis-2rd
### Last motified Dec-8-2022, by Julong wei
### we correct covariates, BrainCV_cell_counts_final_correct_2022-12-06.xlsx
### we directly used the final model decided last time, two datasets analysis together and fit library in the models
### Filtering conditions, at least 30 cells per combination
### autosome genes, CPM>0.5, at least 3 individuals at case and control, respectively
### For each library, at least one case and control individual


outdir <- "./4.0_Diff.cluster.outs/Correct_results/Filter2_default/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


###
MCls_name <- c("0"="ODC", "1"="Astrocyte", "2"="Microglia", "3"="OPC",
     "4"="Non_DA", "6"="Pericyte", "7"="DA", "8"="Endothelial", "9"="T-cell", "12"="Ependymal")
###
MCls_val <- c("ODC"=1, "OPC"=2, "Astrocyte"=3, "Microglia"=4, "DA"=5, "Non_DA"=6, 
    "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)


########################
### pseudo-bulk data ###
########################


### GEO rds
## fn <- "./3_cluster.outs/1.2_seurat.annot.rds"
## sc <- read_rds(fn)
## x <- sc@meta.data
## x2 <- x%>%
##    dplyr::select(orig.ident, nCount_RNA, nFeature_RNA,
##      NEW_BARCODE, DROPLET.TYPE, SNG.BEST.GUESS, EXP, Batch, seurat_clusters, celltype)
## sc <- AddMetaData(sc, x2)

## opfn <- "./3_cluster.outs/SeuratObject-GEO.rds"
## write_rds(sc, file=opfn)



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


###
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


###
### summary results using different threshold #cells per combination
ths <- c(20, 30, 40, 50, 100)

summDF <- map_dfr(ths, function(ii){
###    
ddx <- dd%>%filter(ncell>=ii)
tmp <- ddx%>%mutate(MCls=gsub("_.*", "", bti), sampleID=gsub(".*_", "", bti))%>%as.data.frame()
tmp <- tmp%>%mutate(celltype=MCls_name[as.character(MCls)])
                    
tmp2 <- tmp%>%group_by(celltype)%>%summarise(ncell=n(),.groups="drop")
    
tmp2 <- tmp2%>%mutate(MCls_value=as.numeric(MCls_val[as.character(celltype)]))%>%arrange(MCls_value)

tmp2$ths <- ii
##
tmp2
})

##
summDF2 <- summDF%>%pivot_wider(id_cols=celltype, names_from=ths, names_prefix="ths_", values_from=ncell)

opfn <- "./4.0_Diff.cluster.outs/1_summ.xlsx"
write.xlsx(summDF2, file=opfn, overwrite=T)
    



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

## old covariates 
## x <- read.csv("BrainCV_cell.counts.final.csv")%>%dplyr::filter(USE==1)%>%
##    dplyr::select(sampleID, Library, Age, Race, Sex, BW_lbs, BW_g, pH,
##                  CauseOfDeath, BrainRegion, Category, PC1, PC2, PC3)

###
### correct covariates
x <- openxlsx::read.xlsx("BrainCV_cell_counts_final_correct_2022-12-06.xlsx")%>%
   dplyr::select(sampleID, Library, Age, Race, Sex, BW_lbs, BW_g, pH,
                 CauseOfDeath, BrainRegion, Category, PC1, PC2, PC3)
cvt <- cvt%>%left_join(x, by="sampleID")%>%
    mutate(Data_source=ifelse(grepl("AKB", sampleID), "Bannon", "Mash"))

## cvt_summ <- cvt%>%group_by(MCls, Category)%>%summarise(n=n(),.groups="drop")%>%as.data.frame()%>%
##    mutate(celltype=MCls_name[as.character(MCls)], MCls_value=MCls_val[as.character(celltype)])%>%
##    arrange(MCls_value)%>%           
##    pivot_wider(id_cols=celltype, names_from=Category, values_from=n)

## opfn <- paste(outdir, "2.1_summary.nindi.xlsx", sep="")
## openxlsx::write.xlsx(cvt_summ, file=opfn, overwrite=T)


###
### old covariates
## cvt <- data.frame(bti=bti2, MCls=cvt0[,1], sampleID=cvt0[,2])
## x <- openxlsx::read.xlsx("BrainCV_cell_counts_final_USE_JW.xlsx")%>%
##    dplyr::select(sampleID, Library, Age, Race, Sex, BW_lbs, BW_g, pH,
##                  CauseOfDeath, BrainRegion, Category, PC1, PC2, PC3)
## cvt_old <- cvt%>%left_join(x, by="sampleID")%>%
##     mutate(Data_source=ifelse(grepl("AKB", sampleID), "Bannon", "Mash"))

## cvt_summ <- cvt_old%>%group_by(MCls, Category)%>%summarise(n=n(),.groups="drop")%>%as.data.frame()%>%
##    mutate(celltype=MCls_name[as.character(MCls)], MCls_value=MCls_val[as.character(celltype)])%>%
##    arrange(MCls_value)%>%           
##    pivot_wider(id_cols=celltype, names_from=Category, values_from=n)

## opfn <- paste(outdir, "2.2_old.nindi.xlsx", sep="")
## openxlsx::write.xlsx(cvt_summ, file=opfn, overwrite=T)



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
MCls_name <- c("0"="ODC", "1"="Astrocyte", "2"="Microglia", "3"="OPC", "4"="Non_DA",
   "6"="Pericyte", "7"="DA", "8"="Endothelial", "9"="T-cell")
res <- res%>%mutate(celltype=MCls_name[as.character(MCls)])%>%as.data.frame()

###
### output
opfn <- paste(outdir, "1_DESeq_results_default.rds", sep="")
write_rds(res, opfn)



################
### MA plots ###
################

#########################
### keep all the same ###
#########################

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


###
### read data
## fn <- paste(outdir, "1.2_DESeq.results_CPM0.5_filter.rds", sep="")
## ## fn <- paste(outdir, "1.0_DESeq.results.rds", sep="")
fn <- paste(outdir, "1_DESeq_results_default.rds", sep="")
res <- read_rds(fn)%>%
   mutate(color=ifelse((p.adjusted<0.1)&(!is.na(p.adjusted)), T, F))


figfn <- paste(outdir, "Figure1.1_MA.png", sep="")
png(figfn, width=880, height=500, res=120)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
x <- matrix(1:8, 2, 4, byrow=T)
layout(x)

 
MCls <- c("ODC", "OPC", "Astrocyte", "Microglia", "DA", "Non_DA", "Pericyte")
for (ii in MCls){
   ##1
   ##oneMCl <- MCls[ii] ##paste(ii, MCls_label[ii], sep=":")
   cat(ii, "\n") 
   res2 <- res%>%
       dplyr::filter(celltype==as.character(ii))%>%
       dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
   print(plotMA(res2[,1:3], colLine="NA", main=ii, cex.main=1.2, cex.axis=0.8, cex.lab=1))
}
dev.off()


################
### QQ plots ###
################

fn <- paste(outdir, "1_DESeq_results_default.rds", sep="")
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
    mutate(MCls_value=as.numeric(MCls_val[as.character(celltype)]),
    celltype2=fct_reorder(celltype, MCls_value))
     
### qq plots
p  <- ggplot(plotDF, aes(x=expected, y=observed))+
   geom_point(size=0.3, colour="grey30")+
   facet_wrap(~celltype2, nrow=2, ncol=4, scales="free")+## ,
   ##   labeller=as_labeller(c("0"="0:ODC", "1"="1:astrocyte", "2"="2:microglia", "3"="3:OPC", "4"="4:GABA",
   ## "6"="6:pericyte", "7"="7:DaN", "8"="8:endothelial", "9"="9:T-cells")))+
   geom_abline(colour="red")+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(strip.text=element_text(size=12))
###

figfn <- paste(outdir, "Figure1.2_qq.png", sep="")
png(figfn, width=880, height=500, res=120)
print(p)
dev.off()

###
### qq plots combine cell-types
p2  <- ggplot(plotDF, aes(x=expected, y=observed, col=factor(celltype2)))+
   geom_point(size=0.3)+
   scale_color_manual(values=col_MCls,
       breaks=c("ODC", "OPC", "Astrocyte", "Microglia", "DA", "Non_DA", "Pericyte"),
       guide=guide_legend(override.aes=list(size=2)))+ 
   ## facet_wrap(~celltype2, nrow=2, ncol=4, scales="free")+## ,
   ##   labeller=as_labeller(c("0"="0:ODC", "1"="1:astrocyte", "2"="2:microglia", "3"="3:OPC", "4"="4:GABA",
   ## "6"="6:pericyte", "7"="7:DaN", "8"="8:endothelial", "9"="9:T-cells")))+
   geom_abline(colour="red")+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.key.size=grid::unit(0.8, "lines"))
###

figfn <- paste(outdir, "Figure1.2_qq2_comb.png", sep="")
png(figfn, width=550, height=450, res=120)
print(p2)
dev.off()






#####################
### Table of DEGs ###
#####################
 
outdir <- "./4.0_Diff.cluster.outs/Correct_results/Filter2_default/"    
fn <- paste(outdir, "1_DESeq_results_default.rds", sep="")

## fn <- paste(outdir, "1.2_DESeq.results_CPM0.5_filter.rds", sep="")

res <- read_rds(fn)
res2 <- res%>%mutate(is_DEG=ifelse(p.adjusted<0.1&abs(estimate)>0.25, 1, 0))
 
## opfn <- paste(outdir, "Table1_DEG_list.xlsx", sep="")
## openxlsx::write.xlsx(res2%>%filter(is_DEG==1)%>%dplyr::select(-is_DEG), opfn)

## x1 <- res2%>%dplyr::filter(is_DEG==1, estimate>0)%>%pull(gene)%>%unique()
## x2 <- res2%>%dplyr::filter(is_DEG==1, estimate<0)%>%pull(gene)%>%unique()
 
tmp <- res2%>%group_by(celltype)%>%summarise(ntest=n(), nDEG=sum(is_DEG,na.rm=T), .groups="drop")%>%
    as.data.frame()%>%
    mutate(MCls_value=as.numeric(MCls_val[as.character(celltype)]))%>%
    arrange(MCls_value)

opfn <- paste(outdir, "2.1_summary.DEG.xlsx", sep="")
openxlsx::write.xlsx(tmp, file=opfn, overwrite=T)

############################
### compare Cell reports ###
############################

DEG <- res2%>%filter(celltype=="OPC", is_DEG==1)%>%pull(gene)

x  <- openxlsx::read.xlsx("Cell_reports/Diff_1.xlsx")

DEG_cr <- x$OLs
DEG_cr <- toupper(DEG_cr[!is.na(DEG_cr)]) ##%>%drop_na(Astrocytes)

length(DEG_cr)
ii <- intersect(DEG_cr, DEG)


res2 <- res%>%filter(celltype=="Astrocyte")%>%dplyr::select(gene, x=estimate)

old <- openxlsx::read.xlsx("Cell_reports/Astrocyte_diff.xlsx")%>%
    dplyr::select(gene=GENE, y=log2FC)%>%
    mutate(gene=toupper(gene))


comb <- res2%>%inner_join(old, by="gene")


length(unique(res2$gene)) ## 5239


### old
fn <- "./4.0_Diff.cluster.outs/Old_results/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds"
res <- read_rds(fn)
##
MCls_name <- c("0"="ODC", "1"="Astrocyte", "2"="Microglia", "3"="OPC", "4"="Non_DA",
   "6"="Pericyte", "7"="DA", "8"="Endothelial", "9"="T-cell")
res <- res%>%mutate(celltype=MCls_name[as.character(MCls)])%>%as.data.frame()


res2 <- res%>%mutate(is_DEG=ifelse(p.adjusted<0.1&abs(estimate)>0.25, 1, 0))
tmp <- res2%>%group_by(celltype)%>%summarise(ntest=n(), nDEG=sum(is_DEG,na.rm=T), .groups="drop")%>%
    as.data.frame()%>%
    mutate(MCls_value=as.numeric(MCls_val[as.character(celltype)]))%>%
    arrange(MCls_value)    
 
opfn <- paste(outdir, "2.2_old.DEG.xlsx", sep="")
openxlsx::write.xlsx(tmp, file=opfn, overwrite=T)

res2%>%dplyr::filter(is_DEG==1)%>%pull(gene)%>%unique()%>%length()

length(unique(res2$gene))


#######################################
### summary of differential results ###
#######################################

MCls_val <- c("ODC"=1, "OPC"=2, "Astrocyte"=3, "Microglia"=4, "DA"=5, "Non_DA"=6, 
     "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)

outdir <- "./4.0_Diff.cluster.outs/Correct_results/Filter2_default/"    
fn <- paste(outdir, "1_DESeq_results_default.rds", sep="")

res <- read_rds(fn)
res <- res%>%mutate(is_DEG=ifelse(p.adjusted<0.1&abs(estimate)>0.25, 1, 0),
                    is_UP=ifelse(p.adjusted<0.1&estimate>0.25, 1, 0),
                    is_Down=ifelse(p.adjusted<0.1&estimate<(-0.25), 1, 0))
###
### summary data
summ <- res%>%group_by(celltype)%>%
    summarise(ngene=n(), nDEG=sum(is_DEG, na.rm=T), prop=nDEG/ngene,
        nUP=sum(is_UP, na.rm=T), nDown=sum(is_Down, na.rm=T),
        lambda=median(statistic^2)/qchisq(0.5,1), .groups="drop")%>%ungroup()


###
### old bulk data
fn <- "/wsu/home/groups/bannonlab/manal/DESeq2_WGCNA_DownstreamAnalysis/Final_Runs/DESEQ2_TEST_3-28-18/DESEQ2_IHW-5-28-18.csv"
old <- read.csv(fn)%>%drop_na(pvalue)%>%dplyr::rename("gene"="Associated.Gene.Name")
old_DEG <- old%>%filter(padj<0.1)%>%pull(gene)%>%unique()

res <- res%>%mutate(is_Bulk=ifelse(gene%in%old_DEG, 1, 0))

###
### ks test
MCls <- sort(unique(res$celltype))
resKS <- map_dfr(MCls, function(ii){
   ##
   DF <- res%>%filter(celltype==ii) 
   x1 <- DF%>%filter(is_Bulk==1)%>%pull(p.value)
   x2 <- DF%>%filter(is_Bulk==0)%>%pull(p.value)
   ks <- ks.test(x1, x2, alternative="greater")
   p <- ks$p.value
   score <- ks$statistic
   ##
   resKS <- data.frame(celltype=ii, pval=p, ks_score=score)
   resKS 
})

summ <- summ%>%left_join(resKS, by="celltype")

###
summ <- summ%>%mutate(MCls_value=as.numeric(MCls_val[as.character(celltype)]),
                      celltype2=fct_reorder(celltype, MCls_value))%>%
    arrange(celltype2)

summ2 <- summ[,1:9]

opfn <- paste(outdir, "2.3_summary.DEG.xlsx", sep="")
openxlsx::write.xlsx(summ2, opfn, overwrite=T)

    










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
   dplyr::select(gene=Associated.Gene.Name, zscore_old)     
   ## dplyr::select(gene=Associated.Gene.Name, beta_old=log2FoldChange, zscore_old, pval_old=pvalue) 

### results from sc
## fn <- "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds"
## res2 <- read_rds(fn)%>%drop_na(p.value)%>%
##    mutate(zscore_sc=estimate/stderror)%>%
##    dplyr::select(MCls, gene, beta_sc=estimate, zscore_sc, pval_sc=p.value)

outdir <- "./4.0_Diff.cluster.outs/Correct_results/Filter2_default/"    
fn <- paste(outdir, "1_DESeq_results_default.rds", sep="")
res2 <- read_rds(fn)%>%drop_na(p.value)%>%
    dplyr::select(gene, celltype, zscore_sc=statistic)


plotDF <- res2%>%inner_join(old, by="gene")


###
### scatter plots of zscore between sc-brain and bulk
 
plotDF2 <- plotDF%>%dplyr::select(celltype, gene, x=zscore_sc, y=zscore_old)%>%
    mutate(##celltype=MCls_name[as.character(MCls)],
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
fig1 <- ggplot(plotDF2, aes(x=x, y=y, color=factor(celltype)))+
   geom_point(size=0.3)+
   ##stat_density_2d(aes(fill=..level..), geom="polygon")+ 
   geom_text(data=anno_df1, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)+
   ##scale_fill_viridis(option="C")+ 
   scale_color_manual(values=col_MCls, guide="none")+
   facet_wrap(~celltype2, nrow=2, ncol=4, scales="free")+         
   scale_x_continuous("z-score from different cell-types in sc-brain", expand=expansion(mult=0.1))+
   scale_y_continuous("z-score from Bulk", expand=expansion(mult=0.15))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.position=c(0.9, 0.25),
         legend.key.size=grid::unit(0.8, "lines"),
         strip.text=element_text(size=12))
   
fig1 <- fig1+geom_smooth(method="lm",formula=y~x, size=0.5, se=F, color="red")
                           
figfn <- paste(outdir, "Figure1.3_zscore_SCvsBulk.png", sep="")
png(figfn, width=880, height=500, res=120)
print(fig1)
dev.off()



################################################################################
### scatter plots compare results between correct covariates vs old analysis ###
################################################################################

###
outdir <- "./4.0_Diff.cluster.outs/Correct_results/Filter2_default/"    
fn <- paste(outdir, "1_DESeq_results_default.rds", sep="")
res <- read_rds(fn)%>%drop_na(p.value)%>%
    mutate(comb=paste(celltype, gene, sep="_"),
           is_DEG_correct=ifelse(p.adjusted<0.1&abs(estimate)>0.25, 1, 0))%>%
    dplyr::select(gene, celltype, comb, zscore_correct=statistic, is_DEG_correct)
                             

### 
fn <- "./4.0_Diff.cluster.outs/Old_results/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds"
old <- read_rds(fn)%>%as.data.frame()%>%drop_na(p.value)%>%
    mutate(celltype=MCls_name[as.character(MCls)],
           comb=paste(celltype, gene, sep="_"),
           is_DEG_old=ifelse(p.adjusted<0.1&abs(estimate)>0.25, 1, 0))%>%
    dplyr::select(comb, zscore_old=statistic, is_DEG_old)

plotDF <- res%>%inner_join(old, by="comb")

plotDF <- plotDF%>%
   mutate(gr=case_when( ((is_DEG_correct+is_DEG_old)==0)~"0",
     (is_DEG_correct==1&is_DEG_old==0)~"1",
     (is_DEG_correct==0&is_DEG_old==1)~"2",
     ((is_DEG_correct+is_DEG_old)==2)~"3"))


###
### scatter plots of zscore between correct and old analysis

plotDF2 <- plotDF%>%dplyr::rename(x=zscore_old, y=zscore_correct)%>%
    mutate(MCls_value=as.numeric(MCls_val[as.character(celltype)]),
           celltype2=fct_reorder(celltype, MCls_value))

###
### annotation data frame
anno_df1 <- plotDF2%>%group_by(celltype)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$x, (.x)$y, method="pearson")),
          eq=map(corr,feq),
          r2=map_dbl(corr,~(.x)$estimate),
          xpos=map_dbl(data,~xFun(.x,a=0.4)),
          ypos=map_dbl(data,~yFun(.x,a=0.98)))%>%
   dplyr::select(-data,-corr)%>%unnest(celltype)

anno_df1 <- anno_df1%>%
    mutate(MCls_value=as.numeric(MCls_val[as.character(celltype)]),
           celltype2=fct_reorder(celltype, MCls_value))
 
###
### main figures
p2 <- ggplot(plotDF2, aes(x=x, y=y, color=factor(gr)))+
   geom_point(size=0.3)+
   ##stat_density_2d(aes(fill=..level..), geom="polygon")+ 
   geom_text(data=anno_df1, aes(x=xpos, y=ypos, label=eq), colour="grey30", size=3, parse=T)+
   ##scale_fill_viridis(option="C")+ 
   scale_color_manual(values=c("0"="grey", "1"="red", "2"="blue", "3"="green"),
                      labels=c("0"="not DEG (either)", "1"="DEG (correct)",
                               "2"="DEG (old)", "3"="DEG (both)"),
                      guide=guide_legend(override.aes=list(size=2)))+
   facet_wrap(~celltype2, nrow=2, ncol=4, scales="free")+         
   scale_x_continuous("z-score from old analysis", expand=expansion(mult=0.1))+
   scale_y_continuous("z-score from correct analysis", expand=expansion(mult=0.15))+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.position=c(0.9, 0.3),
         legend.key.size=grid::unit(0.8, "lines"),
         strip.text=element_text(size=12))
   
p2 <- p2+geom_smooth(method="lm",formula=y~x, size=0.5, se=F, color="black")
                           
figfn <- paste(outdir, "Figure1.4_zscore_CorrectVsOld.png", sep="")
png(figfn, width=880, height=500, res=120)
print(p2)
dev.off()





###
###
## fn <-  "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds"
## res <- read_rds(fn)

## summ <- res%>%group_by(MCls)%>%
##    summarise(ntest=sum(!is.na(estimate)),
##              ntest_q=sum(!is.na(p.adjusted)), n_na=sum(is.na(p.adjusted)), .groups="drop")%>%
##    ungroup()

## opfn <- "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/2_summary_qc.xlsx"
## openxlsx::write.xlsx(summ, file=opfn, overwrite=T)


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



   
   
