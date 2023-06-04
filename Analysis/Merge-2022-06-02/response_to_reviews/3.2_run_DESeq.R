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
library(biobroom)
library(qvalue)
##
library(ggrastr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
library(circlize)

rm(list=ls())



outdir <- "./3_re-analysis.outs/Diff_analysis.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


##############################
### ODC run DESeq analysis ###
##############################
 
fn <-  "./3_re-analysis.outs/2.2_Microglia_YtX_sel.comb.rds"
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


###
### correct covariates
x <- openxlsx::read.xlsx("../BrainCV_cell_counts_final_correct_2022-12-06.xlsx")%>%
   dplyr::select(sampleID, Library, Age, Race, Sex, BW_lbs, BW_g, pH,
                 CauseOfDeath, BrainRegion, Category, PC1, PC2, PC3)
cvt <- cvt%>%left_join(x, by="sampleID")%>%
    mutate(Data_source=ifelse(grepl("AKB", sampleID), "Bannon", "Mash"))


###
### run DESeq analysis
th_cpm <- 0.5
MCls <- as.character(sort(as.numeric(unique(cvt$MCls))))
res <- lapply(MCls[1:3], function(ii){

   ### 
   time0 <- Sys.time()
   cvt0 <- cvt%>%dplyr::filter(MCls==ii) 

   ### filter library, at least 1 individual in opioid and control 
   lib_summ <- cvt0%>%group_by(Library, Category)%>%summarise(ny=n(), .groups="drop")%>%
       pivot_wider(id_cols=Library, names_from=Category, values_from=ny, values_fill=0)
   LibSel <- lib_summ%>%dplyr::filter(Control>=1, Opioid>=1)%>%dplyr::pull(Library)
   cvt0 <- cvt0%>%dplyr::filter(Library%in%LibSel)

   
   bti1 <- cvt0%>%filter(Category=="Control")%>%dplyr::pull(bti)
   bti2 <- cvt0%>%filter(Category=="Opioid")%>%dplyr::pull(bti)
   ###
   ### at least 3 individuals in opioid and control across libraries
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
      res2 <- NULL
   }  ### End if try-error

       
   }else{
      res2 <- NULL
   } ### End if >=3
   
   ###
   time1 <- Sys.time()
   elapsed <- difftime(time1, time0, units="mins")
   cat(ii, elapsed, "Done\n")
   
   res2
})

####
###

res <- do.call(rbind,res)
res <- as.data.frame(res)

###
### outputy
opfn <- paste(outdir, "1.2_Microglia_DESeq_results_default.rds", sep="")
write_rds(res, opfn)



################################################
### Summary heatmap of LFC across cell types ###
################################################

###
fn <- "../4.0_Diff.cluster.outs/Correct_results/Filter2_default/1_DESeq_results_default.rds"
res <- read_rds(fn)
resSig <- res%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.25)
DEGs <- unique(resSig$gene)


###
fn <- "./3_re-analysis.outs/Diff_analysis.outs/1.1_ODC_DESeq_results_default.rds"
res2 <- read_rds(fn)
res2 <- res2%>%mutate(celltype=paste("ODC_", MCls, sep=""))


###
fn <- "./3_re-analysis.outs/Diff_analysis.outs/1.2_Microglia_DESeq_results_default.rds"
res3 <- read_rds(fn)
res3 <- res3%>%mutate(celltype=paste("Microglia_", MCls, sep=""))

res3 <- res3%>%group_by(MCls)%>%mutate(FDR=q


###
resAll <- rbind(res, res2, res3)
  
comb <- sort(unique(resAll$celltype))



### correlation matrix for plots
corr <- NULL
icorr <- NULL
for (ii in comb){
    ##
    rr <- NULL
    irr <- NULL
    for (kk in comb){
        ##
        DF <- resAll%>%dplyr::filter(celltype==ii, gene%in%DEGs)%>%dplyr::select(gene, x=estimate)
        DF2 <- resAll%>%dplyr::filter(celltype==kk, gene%in%DEGs)%>%dplyr::select(gene, y=estimate)
        DF <- DF%>%inner_join(DF2, by="gene")%>%drop_na(x,y)
        tmp <- cor.test(DF$x, DF$y, method="pearson")
        rr <- c(rr, tmp$estimate)
        ###
        pval <- tmp$p.value
        is_sig <- ifelse(pval<0.05, 1, 0)
        irr <- c(irr, is_sig)
    } ###
    corr <- cbind(corr, rr)
    icorr <- cbind(icorr, irr)
}
###corr <- corr*icorr
colnames(corr) <- comb
rownames(corr) <- comb

corr <- corr*icorr
###
### colors 
r <- as.numeric(corr)
rneg <- r[r<0]
rpos <- r[r>0]
## mybreak <-c(seq(min(rneg), max(rneg), length.out=10), 0, seq(min(rpos), max(rpos), length.out=10))
mybreak <- seq(-1, 1, length.out=21) 
###
mycol <- colorRamp2(mybreak, colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(21)) 



###
### annotation colors

col_MCls <- c("ODC"="#fb8072", "ODC_0"="#cb181d", "ODC_1"="#fee0d2", "OPC"="#35978f",
   "Astrocyte"="#92cd2d", "Microglia"="#c38dc4", "Microglia_0"="#4a1486", "Microglia_2"="#807dba", 
   "Non_DA"="#fc9016",  "Pericyte"="#bebada", "DA"="#80b1d3", "Endothelial"="#fccde5",
    "T-cell"="#8dd3c7",  "Ependymal"="#828282")

df_col <- data.frame(celltype=colnames(corr))
 
col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col_MCls))

        
###
###
p <- Heatmap(corr, col=mycol, cluster_rows=T, cluster_columns=T,
   row_names_gp=gpar(fontsize=10), column_names_gp=gpar(fontsize=10),
   top_annotation=col_ha,
   heatmap_legend_param=list(title="PCC", legend_height=unit(5,"cm"), grid_width=unit(0.4, "cm")))


figfn <- paste(outdir, "Figure1.0_corr_heatmap.png", sep="")
png(figfn, height=750, width=850, res=120)
set.seed(0)
p <- draw(p)
r.list <- row_order(p)
r.dend <- row_dend(p)
dev.off()


###
### self-defined order

comb_order <- c("ODC", "ODC_0", "ODC_1", "OPC", "Astrocyte",
                "Microglia", "Microglia_0", "Microglia_2",
                "Pericyte", "DA", "Non_DA")
corr2 <- corr[comb_order, comb_order]

### colors
df_col <- data.frame(celltype=colnames(corr2))
col_ha <- HeatmapAnnotation(df=df_col, col=list(celltype=col_MCls))

p2 <- Heatmap(corr2, col=mycol, cluster_rows=F, cluster_columns=F,
   row_names_gp=gpar(fontsize=10), column_names_gp=gpar(fontsize=10),
   top_annotation=col_ha,
   heatmap_legend_param=list(title="PCC", legend_height=unit(5,"cm"), grid_width=unit(0.4, "cm")))


figfn <- paste(outdir, "Figure1.2_corr_heatmap_reorder.png", sep="")
png(figfn, height=650, width=780, res=120)
set.seed(0)
p2 <- draw(p2)
dev.off()


###############################
### summary identified DEGs ###
###############################
###
fn <- "../4.0_Diff.cluster.outs/Correct_results/Filter2_default/1_DESeq_results_default.rds"
res <- read_rds(fn)
resSig <- res%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.25)
DEGs <- unique(resSig$gene)

DEG_ODC <- resSig%>%filter(celltype=="ODC")%>%pull(gene)%>%unique()

###
fn <- "./3_re-analysis.outs/Diff_analysis.outs/1.1_ODC_DESeq_results_default.rds"
res2 <- read_rds(fn)
res2 <- res2%>%mutate(celltype=paste("ODC_", MCls, sep=""))
res2_sig <- res2%>%drop_na(p.value)%>%filter(p.adjusted<0.1, abs(estimate)>0.25)

DEG2_union <- res2_sig%>%pull(gene)%>%unique()
xx <- intersect(DEG_ODC, DEG2_union)

DEG2_cl0 <- res2_sig%>%filter(celltype=="ODC_0")%>%pull(gene)%>%unique()
xx <- intersect(DEG_ODC, DEG2_cl0)

DEG2_cl1 <- res2_sig%>%filter(celltype=="ODC_1")%>%pull(gene)%>%unique()
xx <- intersect(DEG_ODC, DEG2_cl1)

###

DEG_micro <- resSig%>%filter(celltype=="Microglia")%>%pull(gene)%>%unique()

fn <- "./3_re-analysis.outs/Diff_analysis.outs/1.2_Microglia_DESeq_results_default.rds"
res3 <- read_rds(fn)
res3 <- res3%>%mutate(celltype=paste("Microglia_", MCls, sep=""))
res3 <- res3%>%group_by(MCls)%>%mutate(FDR=qvalue(p.value)$qvalues)%>%ungroup()%>%as.data.frame()

res3_sig <- res3%>%drop_na(p.value)%>%filter(p.adjusted<0.1, abs(estimate)>0.25)

 
DEG3_union <- res3_sig%>%pull(gene)%>%unique()
xx <- intersect(DEG_micro, DEG3_union)

DEG3_cl0 <- res3_sig%>%filter(celltype=="Microglia_0")%>%pull(gene)%>%unique()
xx <- intersect(DEG_micro, DEG3_cl0)

DEG3_cl2 <- res3_sig%>%filter(celltype=="Microglia_2")%>%pull(gene)%>%unique()
xx <- intersect(DEG_micro, DEG3_cl2)




