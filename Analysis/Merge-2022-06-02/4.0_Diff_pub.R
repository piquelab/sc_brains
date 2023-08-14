###
library(Matrix)
library(tidyverse)
library(DESeq2)
library(annotables)
library(biobroom)
##library(Seurat)

library(ComplexHeatmap, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(circlize, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(cowplot)
library(ggrepel)
library(ggtext, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(glue)
library(openxlsx)

##

rm(list=ls())

outdir <- "./4.0_Diff.cluster.outs/Correct_results/plots_pub/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)




##############################
#### Plots for publication ###
##############################

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


###################
### UpSet plots ###
###################

### current color
col_MCls <- c("ODC"="#fb8072", "Astrocyte"="#92cd2d", "Microglia"="#c38dc4", "OPC"="#35978f",
   "Non_DA"="#fc9016",  "Pericyte"="#bebada", "DA"="#80b1d3", "Endothelial"="#fccde5",
    "T-cell"="#8dd3c7",  "Ependymal"="#828282")

### old color
## col_MCls <- c("ODC"="#fb8072", "Astrocyte"="#b3de69", "Microglia"="#bebada", "OPC"="#8dd3c7",
##    "Non_DA"="#fdb462",  "Pericyte"="#ffffb3", "DA"="#80b1d3", "Endothelial"="#fccde5",
##     "T-cell"="#35978f", "Ependymal"="#828282")
###
MCls_name <- c("0"="ODC", "1"="Astrocyte", "2"="Microglia", "3"="OPC",
     "4"="Non_DA", "6"="Pericyte", "7"="DA", "8"="Endothelial", "9"="T-cell", "12"="Ependymal")
###
MCls_val <- c("ODC"=1, "OPC"=2, "Astrocyte"=3, "Microglia"=4, "DA"=5, "Non_DA"=6, 
    "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)


## fn <-  "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds"
## res <- read_rds(fn)%>%mutate(celltype=MCls_name[as.character(MCls)])
fn <- "./4.0_Diff.cluster.outs/Correct_results/Filter2_default/1_DESeq_results_default.rds"
res <- read_rds(fn)


### DEGs
res2 <- res%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.25)

###
plot_data <- res2%>%dplyr::select(gene, celltype)
opfn <- paste(outdir, "TableS_Figure2_UpSet.xlsx", sep="")
write.xlsx(plot_data, file=opfn)


celltypes <- sort(unique(res2$celltype))
geneList <- lapply(celltypes, function(ii){
    ##
    DEG <-  res2%>%filter(celltype==ii)%>%pull(gene)%>%unique()
    DEG
})    
names(geneList) <- celltypes


m <- make_comb_mat(geneList, mode="intersect")

m2 <- m[comb_degree(m)>1&comb_size(m)>5]

#0571b0
##
fig <- UpSet(m2,          
   comb_order=order(comb_size(m2), decreasing=T),
   comb_col="#4eb3d3", ##bg_col=c("#7dc7df", "#bce2ee"),
   ##bg_pt_col="#bce2ee", ##pt_size=unit(4, "pt"), lwd=unit(2, "pt"),   
   row_names_gp=gpar(fontsize=10),
   top_annotation=upset_top_annotation(m2, height=unit(4, "cm"), add_numbers=T,
       bar_width=0.8,                                 
       axis_param=list(gp=gpar(fontsize=10)),
       annotation_name_gp=gpar(fontsize=9)),
   right_annotation=upset_right_annotation(m2, add_numbers=T,                                     
     gp=gpar(fill=c("#c38dc4", "#92cd2d", "#35978f", "#fb8072", "#bebada", "#fc9016")),
     axis_param=list(gp=gpar(fontsize=10)),
     annotation_name_gp=gpar(fontsize=10),
     width=unit(5, "cm") )) 
###
figfn <- paste(outdir, "Figure1.2_upset.png", sep="")
png(figfn, width=800, height=580, res=120)
print(fig)
dev.off()


## c("ODC"="#fb8072", "Astrocyte"="#92cd2d", "Microglia"="#c38dc4", "OPC"="#35978f",
##    "Non_DA"="#fc9016",  "Pericyte"="#bebada", "DA"="#80b1d3", "Endothelial"="#fccde5",
##     "T-cell"="#8dd3c7",  "Ependymal"="#828282")

#################################################
### bar plots show overlapping with bulk data ###
#################################################


fn <- "/wsu/home/groups/bannonlab/manal/DESeq2_WGCNA_DownstreamAnalysis/Final_Runs/DESEQ2_TEST_3-28-18/DESEQ2_IHW-5-28-18.csv"
old <- read.csv(fn)%>%drop_na(pvalue)%>%dplyr::rename("gene"="Associated.Gene.Name")
old_DEG <- old%>%filter(padj<0.1)%>%pull(gene)%>%unique()
## length(old_DEG)

celltypes <- sort(unique(res2$celltype))
plotDF <- map_dfr(celltypes, function(ii){
   ###
   gene2 <- res2%>%filter(celltype==ii)%>%pull(gene)%>%unique() 
   ##
   DF <- data.frame(celltype=ii, nolap=length(intersect(gene2, old_DEG)))
   DF
})
plotDF <- rbind(plotDF, data.frame(celltype="DA", nolap=0))

plotDF <- plotDF%>%
    mutate(MCl_value=as.numeric(MCls_val[as.character(celltype)]),
           celltype2=fct_reorder(celltype, MCl_value))


p <- ggplot(plotDF, aes(x=celltype2, y=nolap, fill=celltype))+
   geom_bar(stat="identity", position=position_dodge())+
   geom_text(aes(x=celltype2, y=nolap, label=nolap), vjust=-0.5, size=3)+            
   scale_fill_manual(values=col_MCls, guide="none")+
   ylim(0,200)+  
   ggtitle("Shared DEGs between sc and bulk")+ 
   theme_bw()+
   theme(legend.title=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10),
         axis.text.y=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=12))

figfn <- paste(outdir, "Figure1.3_nolap_barplot.png", sep="")
png(figfn, width=580, height=450, res=120)
print(p)
dev.off()




###
###
celltype <- sort(unique(res2$celltype))
ncell <- length(celltype)
##
plotDF <- NULL
for (i in 1:(ncell-1)){
    ##
    cell1 <- celltype[i]
    res2x <- res2%>%dplyr::filter(celltype==cell1)%>%dplyr::select(gene, beta_x=estimate)
    if ( nrow(res2x)==0) next
    for (j in (i+1):ncell){
        ##
        cell2 <- celltype[j]
        res2y <- res2%>%dplyr::filter(celltype==cell2)%>%dplyr::select(gene, beta_y=estimate)
        ##
        if ( nrow(res2y)==0) next
        
        res2xy <- res2x%>%inner_join(res2y, by="gene")

        if( nrow(res2xy)==0) next
        res2xy <- res2xy%>%mutate(direct=sign(beta_x)*sign(beta_y))
        ##
        DF2 <- res2xy%>%group_by(direct)%>%summarise(ngene=n(),.groups="drop")%>%
            mutate(nt=sum(ngene), comb=paste(cell1, cell2, sep="_"))
        ##
        if ( sum(DF2$ngene)>100){
           plotDF <- rbind(plotDF, DF2)
        }   
    }## cell-2
}## cell-1

###
plotDF <- plotDF%>%mutate(percent=ngene/nt)

 
## #0571b0 
p <- ggplot(plotDF, aes(x=comb, y=percent))+
    geom_bar(stat="identity", aes(alpha=factor(direct)), fill="#0571b0")+   ##"#ca0020")+
    geom_text(aes(label=nt), y=1, vjust=-0.5, size=3)+
    scale_alpha_manual(values=c("-1"=0.4, "1"=1), labels=c("-1"="not same", "1"="same"))+
    ylim(0, 1.1)+ylab("Proportion")+    
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.key.size=unit(0.6, units="cm"),
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=-45, hjust=0, vjust=0.5))
         # plot.margin=grid::unit(c(5.5,5.5, 5.5, 6), "points"))
          ##axis.text.y=element_text(size=10))

###
figfn <- paste(outdir, "Figure1.4_olap_direction.png", sep="")
png(figfn, width=580, height=420, res=120)
p
dev.off()






################
### QQ plots ###
################
 
myks <-function(DF){
   ##
   x1 <- DF%>%filter(is_bulk==1)%>%pull(p.value)
   x2 <- DF%>%filter(is_bulk==2)%>%pull(p.value)
   ks <- ks.test(x1, x2, alternative="greater")
   p <- ks$p.value
   score <- round(as.numeric(ks$statistic),digits=3)

   ##
   if(p<0.001) symb <- "***"
   if(p>=0.001 & p<0.01) symb <- "**"
   if (p>=0.01 & p<0.05) symb <- "*"
   if(p>0.05) symb <- "NS"
   pval <- round(as.numeric(p), digits=3) 
   ### 
   eq <- bquote("K-S test"~.(symb))
   eq
}    


##
xFun <- function(dx,a=0.5){
min1 <- min(dx$expected, na.rm=T)
max2 <- max(dx$expected, na.rm=T)
R <- max2-min1
xpos <- min1+a*R
}
##
yFun <- function(dx, a=0.8){
min1 <- min(dx$observed, na.rm=T)
max2 <- max(dx$observed, na.rm=T)
R <- max2-min1
ypos <- min1+a*R
}


### setting parameters
MCls_name <- c("0"="ODC", "1"="Astrocyte", "2"="Microglia", "3"="OPC",
     "4"="Non_DA", "6"="Pericyte", "7"="DA", "8"="Endothelial", "9"="T-cell", "12"="Ependymal")

## old color 
## col_MCls <- c("ODC"="#fb8072", "Astrocyte"="#b3de69", "Microglia"="#bebada", "OPC"="#8dd3c7",
##    "Non_DA"="#fdb462",  "Pericyte"="#ffffb3", "DA"="#80b1d3", "Endothelial"="#fccde5",
##     "T-cell"="#35978f", "Ependymal"="#828282")

MCls_val <- c("ODC"=1, "OPC"=2, "Astrocyte"=3, "Microglia"=4, "DA"=5, "Non_DA"=6, 
    "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)




###
### plot data

## fn <-  "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds"
## res <- read_rds(fn)%>%drop_na(p.value)%>%as.data.frame()
fn <- "./4.0_Diff.cluster.outs/Correct_results/Filter2_default/1_DESeq_results_default.rds"
res <- read_rds(fn)%>%drop_na(p.value)



###
### old bulk data
fn <- "/wsu/home/groups/bannonlab/manal/DESeq2_WGCNA_DownstreamAnalysis/Final_Runs/DESEQ2_TEST_3-28-18/DESEQ2_IHW-5-28-18.csv"
old <- read.csv(fn)%>%drop_na(pvalue)%>%dplyr::rename("gene"="Associated.Gene.Name")
old_DEG <- old%>%filter(padj<0.1)%>%pull(gene)%>%unique()

 
res <- res%>%mutate(is_bulk=ifelse(gene%in%old_DEG, 1, 2))
res <- res%>%mutate(celltype=as.character(MCls_name[as.character(MCls)]))

res2 <- res%>%group_by(celltype, is_bulk)%>%arrange(p.value)%>%
   mutate(observed=-log10(p.value), expected=-log10(ppoints(length(p.value))))%>%
   ungroup()%>%as.data.frame() 


plotDF <- res2%>%
   mutate(MCls_value=as.numeric(MCls_val[as.character(celltype)]),
          celltype2=fct_reorder(celltype, as.numeric(MCls_value)))



anno_df <- plotDF%>%group_by(celltype)%>%
   nest()%>%
   mutate(eq=map(data, myks), xpos=1.5, ypos=11)%>%
   dplyr::select(-data)

anno_df <- anno_df%>%
    mutate(MCls_value=as.numeric(MCls_val[celltype]),
           celltype2=fct_reorder(celltype, MCls_value))

        
###
###

p <- ggplot(plotDF)+
   geom_point(aes(x=expected, y=observed, colour=factor(is_bulk)), size=0.3)+
   geom_text(data=anno_df, aes(x=xpos, y=ypos, label=eq), size=3, parse=T)+ 
   scale_colour_manual(values=c("1"="#e7298a", "2"="#1b9e77"),
       labels=c("1"="DEGs in Bulk", "2"="Not DEGs"),
       guide=guide_legend(override.aes=list(size=2)))+
   ## scale_colour_manual(values=col_MCls,
   ##     breaks=c("ODC", "Astrocyte", "Microglia", "OPC", "Non_DA", "Pericyte", "DA"),
   ##     guide=guide_legend(override.aes=list(size=2)))+
   ## scale_alpha_manual(values=c("1"=1, "2"=0.5), labels=c("1"="DEGs in Bulk", "2"="Not DEGs"),
   ##     guide=guide_legend(override.aes=list(size=2, colour="#b3de69")))+
   geom_abline(intercept=0, slope=1)+
   xlab(bquote("Expected"~-log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+ylim(0,13)+
   facet_wrap(~celltype2, nrow=2, ncol=4, scales="fixed")+ 
   ## facet_wrap(~, labeller=as_labeller(c("1"="DEG in Bulk", "2"="Not DEG in Bulk")), scales="free", nrow=2)+ 
   theme_bw()+
   theme(legend.title=element_blank(),
         #legend.text=element_text(size=10),
         legend.key.size=grid::unit(0.8, "lines"),
         legend.position=c(0.9,0.25),
         strip.text=element_text(size=12))
         #axis.text=element_text(size=10),
         #axis.title=element_text(size=10))
 
figfn <- paste(outdir, "Figure2.1_DEG_qq.png", sep="")
png(figfn, width=880, height=500, res=120)
print(p)
dev.off()




####################################################################
### barplots show the correlation of LFC with previous bulk data ###
####################################################################

### function
feq <- function(x){
  r <- round(as.numeric(x$estimate),digits=3)
  p <- x$p.value
  if(p<0.001) symb <- "***"
  if(p>=0.001 & p<0.01) symb <- "**"
  if (p>=0.01 & p<0.05) symb <- "*"
  if(p>0.05) symb <- "NS"
  
  eq <- bquote(.(r)~.(symb))
  eq 
}



## current color
col_MCls <- c("ODC"="#fb8072", "Astrocyte"="#92cd2d", "Microglia"="#c38dc4", "OPC"="#35978f",
   "Non_DA"="#fc9016",  "Pericyte"="#bebada", "DA"="#80b1d3", "Endothelial"="#fccde5",
    "T-cell"="#8dd3c7",  "Ependymal"="#828282")
### setting parameters
MCls_name <- c("0"="ODC", "1"="Astrocyte", "2"="Microglia", "3"="OPC",
     "4"="Non_DA", "6"="Pericyte", "7"="DA", "8"="Endothelial", "9"="T-cell", "12"="Ependymal")


MCls_val <- c("ODC"=1, "OPC"=2, "Astrocyte"=3, "Microglia"=4, "DA"=5, "Non_DA"=6, 
    "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)


### sc results
fn <- "./4.0_Diff.cluster.outs/Correct_results/Filter2_default/1_DESeq_results_default.rds"
res <- read_rds(fn)%>%drop_na(p.value)

## fn <-  "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds"
## res <- read_rds(fn)%>%drop_na(p.value)%>%
##    mutate(celltype=MCls_name[as.character(MCls)])%>%as.data.frame()
###
res2 <- res%>%dplyr::select(gene, zscore_sc=statistic, MCls, celltype)

### DEGs
##DEG_sc <- res%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.25)%>%pull(gene)%>%unique()



###
### old bulk data
fn <- "/wsu/home/groups/bannonlab/manal/DESeq2_WGCNA_DownstreamAnalysis/Final_Runs/DESEQ2_TEST_3-28-18/DESEQ2_IHW-5-28-18.csv"
old <- read.csv(fn)%>%drop_na(pvalue)%>%dplyr::rename("gene"="Associated.Gene.Name")
old2 <- old%>%drop_na(stat)%>%dplyr::select(gene, zscore_bulk=stat)
##DEG_old <- old%>%filter(padj<0.1)%>%pull(gene)%>%unique()
## length(old_DEG)

##
comb <- res2%>%inner_join(old2, by="gene")

##
plotDF <- comb%>%group_by(celltype)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$zscore_sc, (.x)$zscore_bulk, method="pearson")),
          eq=map(corr,feq), r2=map_dbl(corr,~(.x)$estimate), pval=map_dbl(corr, ~(.x)$p.value))%>%
   dplyr::select(-data,-corr)%>%as.data.frame()
 
plotDF <- plotDF%>%
    mutate(MCl_value=as.numeric(MCls_val[as.character(celltype)]),
           celltype2=fct_reorder(celltype, MCl_value))


p <- ggplot(plotDF, aes(x=celltype2, y=r2, fill=celltype))+
   geom_bar(stat="identity", position=position_dodge())+
   geom_text(aes(x=celltype2, y=r2, label=eq), vjust=-0.5, size=3, parse=T)+            
   scale_fill_manual(values=col_MCls, guide="none")+
   ylim(0, 0.35)+
   ggtitle("PCC of standardized LFC (sc and bulk)")+ 
   theme_bw()+
   theme(legend.title=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         axis.text.x=element_text(angle=45, hjust=1),
         ##axis.text.y=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=12))

figfn <- paste(outdir, "Figure3.1_r2_barplot.png", sep="")
png(figfn, width=580, height=450, res=120)
print(p)
dev.off()



### significant DEG in bulk
### plot data
## plotDF <- comb%>%filter(gene%in%DEG_old)%>%group_by(celltype)%>%
##    nest()%>%
##    mutate(corr=map(data, ~cor.test((.x)$zscore_sc, (.x)$zscore_bulk, method="pearson")),
##           eq=map(corr,feq), r2=map_dbl(corr,~(.x)$estimate))%>%
##    dplyr::select(-data,-corr)%>%as.data.frame()

## plotDF <- plotDF%>%
##     mutate(MCl_value=as.numeric(MCls_val[as.character(celltype)]),
##            celltype2=fct_reorder(celltype, MCl_value))


## p2 <- ggplot(plotDF, aes(x=celltype2, y=r2, fill=celltype))+
##    geom_bar(stat="identity", position=position_dodge())+
##    geom_text(aes(x=celltype2, y=r2, label=eq), vjust=-0.5, size=3, parse=T)+            
##    scale_fill_manual(values=col_MCls, guide="none")+
##    ylim(0, 0.55)+
##    ggtitle("Correlation of standardized LFC (sc and bulk)")+  
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          axis.title.x=element_blank(),
##          axis.title=element_blank(),
##          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10),
##          axis.text.y=element_text(size=10),
##          plot.title=element_text(hjust=0.5, size=12))

## figfn <- paste(outdir, "Figure3.2_r2_barplot.png", sep="")
## png(figfn, width=580, height=450, res=120)
## print(p2)
## dev.off()   

   

##################################
### Forest plots show examples ###
##################################


## current color
col_MCls <- c("ODC"="#fb8072", "Astrocyte"="#92cd2d", "Microglia"="#c38dc4", "OPC"="#35978f",
   "Non_DA"="#fc9016",  "Pericyte"="#bebada", "DA"="#80b1d3", "Endothelial"="#fccde5",
    "T-cell"="#8dd3c7",  "Ependymal"="#828282", "Bulk"="grey40")
### setting parameters
MCls_name <- c("0"="ODC", "1"="Astrocyte", "2"="Microglia", "3"="OPC",
     "4"="Non_DA", "6"="Pericyte", "7"="DA", "8"="Endothelial", "9"="T-cell", "12"="Ependymal")


MCls_val <- c("ODC"=1, "OPC"=2, "Astrocyte"=3, "Microglia"=4, "DA"=5, "Non_DA"=6, 
    "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10, "Bulk"=11)




## Mycol2 <- c(col_MCls, "Bulk"="grey40", c("0"="#fb8072", "1"="#b3de69", "2"="#bebada", "3"="#8dd3c7",
##    "4"="#fdb462",  "6"="#ffffb3", "7"="#80b1d3", "8_Bulk"="#ca0020",
##    "not_DEG"="grey40")
 
## fn <-  "./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds"
## res <- read_rds(fn)%>%drop_na(p.value)%>%as.data.frame()%>%
##    mutate(celltype=MCls_name[as.character(MCls)])


fn <- "./4.0_Diff.cluster.outs/Correct_results/Filter2_default/1_DESeq_results_default.rds"
res <- read_rds(fn)%>%drop_na(p.value)


###
### old bulk data
fn <- "/wsu/home/groups/bannonlab/manal/DESeq2_WGCNA_DownstreamAnalysis/Final_Runs/DESEQ2_TEST_3-28-18/DESEQ2_IHW-5-28-18.csv"
old <- read.csv(fn)%>%drop_na(pvalue)%>%dplyr::rename("gene"="Associated.Gene.Name")

## x2 <- old%>%filter(padj<0.1)%>%arrange(desc(abs(log2FoldChange)))%>%
##     dplyr::select(baseMean, log2FoldChange, lfcSE, gene, pvalue, padj)
## length(old_DEG)


##
## motify data frame

my_plotDF <- function(plotDF){
###
   genes <- sort(unique(plotDF$gene))

   ### 
   plotDF2 <- NULL
   celltypes <- sort(unique(plotDF$celltype))
   for (ii in genes){
       ##
       DF <- plotDF%>%filter(gene==ii)
       ###
       ### 
       celltype2 <- setdiff(celltypes, unique(DF$celltype))
       ncell2 <- length(celltype2) 
       if ( ncell2>0){
          ## 
          DF2 <- data.frame(gene=rep(ii, ncell2), beta=rep(NA, ncell2), SE=rep(NA, ncell2),
                            p.adjusted=rep(NA, ncell2), celltype=celltype2)
          DF <- rbind(DF, DF2)
       }
       plotDF2 <- rbind(plotDF2, DF)
    }
    plotDF2
}    
       


#############################################################################################
### 16 genes including  10 previous DEGs reproted in paper and 6 new genes Mike suggested ###
#############################################################################################


 
DEGlist <- c("PLA1A", "MAP3K6", "YBX3", "FOSL2", "MGP", "MEFV", "CISH", "PDLIM1", "GPR4", "TRIP10",
    "IL4R", "CDKN1A", "CCL2", "NFKBIA", "MIR210HG", "A2M")
###
res2 <- res%>%filter(gene%in%DEGlist)%>%
   dplyr::select(gene, beta=estimate, SE=stderror, p.adjusted, celltype)

###
old2 <- old%>%filter(gene%in%DEGlist)%>%mutate(celltype="Bulk")%>%
   dplyr::select(gene, beta=log2FoldChange, SE=lfcSE, p.adjusted=padj, celltype)



###
plotDF <- rbind(res2, old2)
plotDF2 <- my_plotDF(plotDF)

plotDF2 <- plotDF2%>%
   mutate(p.adjusted=ifelse(is.na(p.adjusted), 1, p.adjusted),
          celltypeNew=ifelse(p.adjusted<0.1&abs(beta)>0, glue("<b>{celltype}"), glue("<i>{celltype}")),
          MCls_value=as.numeric(MCls_val[as.character(celltype)]),
          celltypeNew2=fct_reorder(celltypeNew, MCls_value),
          beta_upper=beta+1.96*SE,
          beta_lower=beta-1.96*SE)


###
### main
p0 <- ggplot(plotDF2, aes(x=beta, y=factor(celltypeNew2), color=factor(celltype)))+
   geom_point(shape=19, size=1.5, na.rm=T)+    
   geom_errorbarh(aes(xmax=beta_upper, xmin=beta_lower), size=0.5, height=0.2, na.rm=T)+
   scale_colour_manual(values=col_MCls, guide="none")+
   xlab("Log2FoldChange")+
   ##scale_y_discrete(breaks=c("ODC", "OPC", "Astrocyte", "Microglia", "DA", "Non_DA", "Pericyte", "Bulk"))+ 
   geom_vline(xintercept=0, size=0.25, linetype="dashed")+
   facet_wrap(~factor(gene), ncol=4, scales="free")+
   theme_bw()+
   theme(strip.text=element_text(size=12),
         axis.title.y=element_blank(),
         axis.text.y=element_markdown())


figfn <- paste(outdir, "Figure4.0_forestplot.png", sep="")
png(filename=figfn, width=1000, height=1000, res=120)  
print(p0)
dev.off()



###
### source data

plotDF2 <- plotDF2%>%
   mutate(p.adjusted=ifelse(is.na(p.adjusted), 1, p.adjusted),
          celltypeNew=gsub("_", " ", celltype),
          ##celltypeNew=ifelse(p.adjusted<0.1&abs(beta)>0, glue("<b>{celltype}"), glue("<i>{celltype}")),
          MCls_value=as.numeric(MCls_val[as.character(celltype)]),
          ##celltypeNew2=fct_reorder(celltypeNew, MCls_value),
          beta_upper=beta+1.96*SE,
          beta_lower=beta-1.96*SE)

opfn <- paste(outdir, "TableS_Figure4_forestplot.txt", sep="")
write.table(plotDF2, opfn, row.names=F, col.names=T, quote=F, sep="\t")



#########################################################################
### select top 10 genes by p-value in sc but not significant in bulk  ###
#########################################################################


DEG_sc <- res%>%dplyr::filter(p.adjusted<0.1, abs(estimate)>0.25)%>%dplyr::pull(gene)%>%unique()
DEG_old <- old%>%dplyr::filter(padj<0.1)%>%pull(gene)%>%unique()
DEG_unq <- setdiff(DEG_sc, DEG_old)


### top 10 DEGs
DEG_top <- res%>%dplyr::filter(gene%in%DEG_unq)%>%
    group_by(gene)%>%
    slice_max(order_by=abs(statistic), n=1)%>%ungroup()%>%
    arrange(desc(abs(statistic)))%>%pull(gene)%>%unique()
DEGlist <- DEG_top[1:10]


###
res2 <- res%>%filter(gene%in%DEGlist)%>%
   dplyr::select(gene, beta=estimate, SE=stderror, p.adjusted, celltype)

###
old2 <- old%>%filter(gene%in%DEGlist)%>%mutate(celltype="Bulk")%>%
   dplyr::select(gene, beta=log2FoldChange, SE=lfcSE, p.adjusted=padj, celltype)



###
plotDF <- rbind(res2, old2)
plotDF2 <- my_plotDF(plotDF)

plotDF2 <- plotDF2%>%
   mutate(p.adjusted=ifelse(is.na(p.adjusted), 1, p.adjusted),
          celltypeNew=ifelse(p.adjusted<0.1&abs(beta)>0, glue("<b>{celltype}"), glue("<i>{celltype}")),
          MCls_value=as.numeric(MCls_val[as.character(celltype)]),
          celltypeNew2=fct_reorder(celltypeNew, MCls_value),
          beta_upper=beta+1.96*SE,
          beta_lower=beta-1.96*SE)


###
### main
p1 <- ggplot(plotDF2, aes(x=beta, y=factor(celltypeNew2), color=factor(celltype)))+
   geom_point(shape=19, size=1.5, na.rm=T)+    
   geom_errorbarh(aes(xmax=beta_upper, xmin=beta_lower), size=0.5, height=0.2, na.rm=T)+
   scale_colour_manual(values=col_MCls, guide="none")+
   xlab("Log2FoldChange")+
   ##scale_y_discrete(breaks=c("ODC", "OPC", "Astrocyte", "Microglia", "DA", "Non_DA", "Pericyte", "Bulk"))+ 
   geom_vline(xintercept=0, size=0.25, linetype="dashed")+
   facet_wrap(~factor(gene), ncol=5, scales="free")+
   theme_bw()+
   theme(strip.text=element_text(size=12),
         axis.title.y=element_blank(),
         ##axis.text.x=element_text(size=8),
         axis.text.y=element_markdown(),
         plot.margin=unit(c(0.1, 0.5, 0.1, 0.1), "cm"))

  
figfn <- paste(outdir, "Figure4.1_forestplot.png", sep="")
png(filename=figfn, width=1100, height=500, res=120)  
print(p1)
dev.off()




###
### TWAS

fn <- "/wsu/home/groups/bannonlab/sc_brains/Analysis/TWAS/1_TWAS.outs/1_TWAS_traits.xlsx"
df <- openxlsx::read.xlsx(fn)
fn <- "/wsu/home/groups/bannonlab/sc_brains/Analysis/TWAS/1_TWAS.outs/traits.xlsx"
trait_DF <- openxlsx::read.xlsx(fn)%>%dplyr::select(traits, Relevant)
df <- df%>%left_join(trait_DF, by="traits")%>%filter(Relevant==1)
twas_gene <- unique(df$SYMBOL)

DEG_unq <- intersect(setdiff(DEG_sc, DEG_old), twas_gene)


DEG_top <- res%>%dplyr::filter(gene%in%DEG_unq)%>%
    group_by(gene)%>%
    slice_max(order_by=abs(statistic), n=1)%>%ungroup()%>%
    arrange(desc(abs(statistic)))%>%pull(gene)%>%unique()
DEGlist <- DEG_top[1:10]


###
res2 <- res%>%filter(gene%in%DEGlist)%>%
   dplyr::select(gene, beta=estimate, SE=stderror, p.adjusted, celltype)

###
old2 <- old%>%filter(gene%in%DEGlist)%>%mutate(celltype="Bulk")%>%
   dplyr::select(gene, beta=log2FoldChange, SE=lfcSE, p.adjusted=padj, celltype)



###
plotDF <- rbind(res2, old2)
plotDF2 <- my_plotDF(plotDF)

plotDF2 <- plotDF2%>%
   mutate(p.adjusted=ifelse(is.na(p.adjusted), 1, p.adjusted),
          celltypeNew=ifelse(p.adjusted<0.1&abs(beta)>0, glue("<b>{celltype}"), glue("<i>{celltype}")),
          MCls_value=as.numeric(MCls_val[as.character(celltype)]),
          celltypeNew2=fct_reorder(celltypeNew, MCls_value),
          beta_upper=beta+1.96*SE,
          beta_lower=beta-1.96*SE)


###
### main
p2 <- ggplot(plotDF2, aes(x=beta, y=factor(celltypeNew2), color=factor(celltype)))+
   geom_point(shape=19, size=1.5, na.rm=T)+    
   geom_errorbarh(aes(xmax=beta_upper, xmin=beta_lower), size=0.5, height=0.2, na.rm=T)+
   scale_colour_manual(values=col_MCls, guide="none")+
   xlab("Log2FoldChange")+
   ##scale_y_discrete(breaks=c("ODC", "OPC", "Astrocyte", "Microglia", "DA", "Non_DA", "Pericyte", "Bulk"))+ 
   geom_vline(xintercept=0, size=0.25, linetype="dashed")+
   facet_wrap(~factor(gene), ncol=5, scales="free")+
   theme_bw()+
   theme(strip.text=element_text(size=12),
         axis.title.y=element_blank(),
         ##axis.text.x=element_text(size=8),
         axis.text.y=element_markdown(),
         plot.margin=unit(c(0.1, 0.5, 0.1, 0.1), "cm"))
  
figfn <- paste(outdir, "Figure4.2_forestplot.png", sep="")
png(filename=figfn, width=1100, height=500, res=120)  
print(p2)
dev.off()


### 
## x <- df%>%filter(SYMBOL%in%DEGlist)%>%arrange(SYMBOL)






#####################################
#### TWAS genes for each category ###
#####################################

fn <- "/wsu/home/groups/bannonlab/sc_brains/Analysis/TWAS/1_TWAS.outs/1_TWAS_traits.xlsx"
df <- openxlsx::read.xlsx(fn)
fn <- "/wsu/home/groups/bannonlab/sc_brains/Analysis/TWAS/1_TWAS.outs/traits.xlsx"
trait_DF <- openxlsx::read.xlsx(fn)%>%dplyr::select(traits, Relevant)
df <- df%>%left_join(trait_DF, by="traits")%>%filter(Relevant==1)

##
## gene list, significant in SC but not in bulk and also related to traits of interest
## DEG_unq


###
### category ADHD
ii <- "ADHD"
geneSel <- df%>%filter(SYMBOL%in%DEG_unq, category==ii)%>%pull(SYMBOL)%>%unique()

###
res2 <- res%>%filter(gene%in%geneSel)%>%
   dplyr::select(gene, beta=estimate, SE=stderror, p.adjusted, celltype)

###
old2 <- old%>%filter(gene%in%geneSel)%>%mutate(celltype="Bulk")%>%
   dplyr::select(gene, beta=log2FoldChange, SE=lfcSE, p.adjusted=padj, celltype)



###
plotDF <- rbind(res2, old2)
plotDF2 <- my_plotDF(plotDF)

plotDF2 <- plotDF2%>%
   mutate(p.adjusted=ifelse(is.na(p.adjusted), 1, p.adjusted),
          celltypeNew=ifelse(p.adjusted<0.1&abs(beta)>0, glue("<b>{celltype}"), glue("<i>{celltype}")),
          MCls_value=as.numeric(MCls_val[as.character(celltype)]),
          celltypeNew2=fct_reorder(celltypeNew, MCls_value),
          beta_upper=beta+1.96*SE,
          beta_lower=beta-1.96*SE)


###
### main
p3 <- ggplot(plotDF2, aes(x=beta, y=factor(celltypeNew2), color=factor(celltype)))+
   geom_point(shape=19, size=1.5, na.rm=T)+    
   geom_errorbarh(aes(xmax=beta_upper, xmin=beta_lower), size=0.5, height=0.2, na.rm=T)+
   scale_colour_manual(values=col_MCls, guide="none")+
   xlab("Log2FoldChange")+
   ##scale_y_discrete(breaks=c("ODC", "OPC", "Astrocyte", "Microglia", "DA", "Non_DA", "Pericyte", "Bulk"))+ 
   geom_vline(xintercept=0, size=0.25, linetype="dashed")+
   facet_wrap(~factor(gene), ncol=3, scales="free")+
   theme_bw()+
   theme(strip.text=element_text(size=12),
         axis.title.y=element_blank(),
         ##axis.text.x=element_text(size=8),
         axis.text.y=element_markdown())
         ##plot.margin=unit(c(0.1, 0.5, 0.1, 0.1), "cm"))
  
figfn <- paste(outdir, "Figure4.3_", ii, "_forestplot.png", sep="")
png(filename=figfn, width=600, height=280, res=120)  
print(p3)
dev.off()


###
### parkinson

ii <- "parkinson"
geneSel <- df%>%filter(SYMBOL%in%DEG_unq, category==ii)%>%pull(SYMBOL)%>%unique()

###
res2 <- res%>%filter(gene%in%geneSel)%>%
   dplyr::select(gene, beta=estimate, SE=stderror, p.adjusted, celltype)

###
old2 <- old%>%filter(gene%in%geneSel)%>%mutate(celltype="Bulk")%>%
   dplyr::select(gene, beta=log2FoldChange, SE=lfcSE, p.adjusted=padj, celltype)



###
plotDF <- rbind(res2, old2)
plotDF2 <- my_plotDF(plotDF)

plotDF2 <- plotDF2%>%
   mutate(p.adjusted=ifelse(is.na(p.adjusted), 1, p.adjusted),
          celltypeNew=ifelse(p.adjusted<0.1&abs(beta)>0, glue("<b>{celltype}"), glue("<i>{celltype}")),
          MCls_value=as.numeric(MCls_val[as.character(celltype)]),
          celltypeNew2=fct_reorder(celltypeNew, MCls_value),
          beta_upper=beta+1.96*SE,
          beta_lower=beta-1.96*SE)


###
### main
p3 <- ggplot(plotDF2, aes(x=beta, y=factor(celltypeNew2), color=factor(celltype)))+
   geom_point(shape=19, size=1.5, na.rm=T)+    
   geom_errorbarh(aes(xmax=beta_upper, xmin=beta_lower), size=0.5, height=0.2, na.rm=T)+
   scale_colour_manual(values=col_MCls, guide="none")+
   xlab("Log2FoldChange")+
   ##scale_y_discrete(breaks=c("ODC", "OPC", "Astrocyte", "Microglia", "DA", "Non_DA", "Pericyte", "Bulk"))+ 
   geom_vline(xintercept=0, size=0.25, linetype="dashed")+
   facet_wrap(~factor(gene), ncol=2, scales="free")+
   theme_bw()+
   theme(strip.text=element_text(size=12),
         axis.title.y=element_blank(),
         ##axis.text.x=element_text(size=8),
         axis.text.y=element_markdown())
         ##plot.margin=unit(c(0.1, 0.5, 0.1, 0.1), "cm"))
  
figfn <- paste(outdir, "Figure4.3_", ii, "_forestplot.png", sep="")
png(filename=figfn, width=480, height=300, res=120)  
print(p3)
dev.off()


###
### marijuana

ii <- "marijuana"
geneSel <- df%>%filter(SYMBOL%in%DEG_unq, category==ii)%>%pull(SYMBOL)%>%unique()

###
res2 <- res%>%filter(gene%in%geneSel)%>%
   dplyr::select(gene, beta=estimate, SE=stderror, p.adjusted, celltype)

###
old2 <- old%>%filter(gene%in%geneSel)%>%mutate(celltype="Bulk")%>%
   dplyr::select(gene, beta=log2FoldChange, SE=lfcSE, p.adjusted=padj, celltype)



###
plotDF <- rbind(res2, old2)
plotDF2 <- my_plotDF(plotDF)

plotDF2 <- plotDF2%>%
   mutate(p.adjusted=ifelse(is.na(p.adjusted), 1, p.adjusted),
          celltypeNew=ifelse(p.adjusted<0.1&abs(beta)>0, glue("<b>{celltype}"), glue("<i>{celltype}")),
          MCls_value=as.numeric(MCls_val[as.character(celltype)]),
          celltypeNew2=fct_reorder(celltypeNew, MCls_value),
          beta_upper=beta+1.96*SE,
          beta_lower=beta-1.96*SE)


###
### main
p3 <- ggplot(plotDF2, aes(x=beta, y=factor(celltypeNew2), color=factor(celltype)))+
   geom_point(shape=19, size=1.5, na.rm=T)+    
   geom_errorbarh(aes(xmax=beta_upper, xmin=beta_lower), size=0.5, height=0.2, na.rm=T)+
   scale_colour_manual(values=col_MCls, guide="none")+
   xlab("Log2FoldChange")+
   ##scale_y_discrete(breaks=c("ODC", "OPC", "Astrocyte", "Microglia", "DA", "Non_DA", "Pericyte", "Bulk"))+ 
   geom_vline(xintercept=0, size=0.25, linetype="dashed")+
   facet_wrap(~factor(gene), ncol=4, scales="free")+
   theme_bw()+
   theme(strip.text=element_text(size=12),
         axis.title.y=element_blank(),
         ##axis.text.x=element_text(size=8),
         axis.text.y=element_markdown())
         ##plot.margin=unit(c(0.1, 0.5, 0.1, 0.1), "cm"))
  
figfn <- paste(outdir, "Figure4.3_", ii, "_forestplot.png", sep="")
png(filename=figfn, width=900, height=500, res=120)  
print(p3)
dev.off()


###
### caffeine
 
ii <- "caffeine"
geneSel <- df%>%filter(SYMBOL%in%DEG_unq, category==ii)%>%pull(SYMBOL)%>%unique()

DEG_top <- res%>%dplyr::filter(gene%in%geneSel)%>%
    group_by(gene)%>%
    slice_max(order_by=abs(statistic), n=1)%>%ungroup()%>%
    arrange(desc(abs(statistic)))%>%pull(gene)%>%unique()


###
res2 <- res%>%filter(gene%in%as.character(DEG_top[1:10]))%>%
   dplyr::select(gene, beta=estimate, SE=stderror, p.adjusted, celltype)

###
old2 <- old%>%filter(gene%in%as.character(DEG_top[1:10]))%>%mutate(celltype="Bulk")%>%
   dplyr::select(gene, beta=log2FoldChange, SE=lfcSE, p.adjusted=padj, celltype)



###
plotDF <- rbind(res2, old2)
plotDF2 <- my_plotDF(plotDF)

plotDF2 <- plotDF2%>%
   mutate(p.adjusted=ifelse(is.na(p.adjusted), 1, p.adjusted),
          celltypeNew=ifelse(p.adjusted<0.1&abs(beta)>0, glue("<b>{celltype}"), glue("<i>{celltype}")),
          MCls_value=as.numeric(MCls_val[as.character(celltype)]),
          celltypeNew2=fct_reorder(celltypeNew, MCls_value),
          beta_upper=beta+1.96*SE,
          beta_lower=beta-1.96*SE)


###
### main
p3 <- ggplot(plotDF2, aes(x=beta, y=factor(celltypeNew2), color=factor(celltype)))+
   geom_point(shape=19, size=1.5, na.rm=T)+    
   geom_errorbarh(aes(xmax=beta_upper, xmin=beta_lower), size=0.5, height=0.2, na.rm=T)+
   scale_colour_manual(values=col_MCls, guide="none")+
   xlab("Log2FoldChange")+
   ##scale_y_discrete(breaks=c("ODC", "OPC", "Astrocyte", "Microglia", "DA", "Non_DA", "Pericyte", "Bulk"))+ 
   geom_vline(xintercept=0, size=0.25, linetype="dashed")+
   facet_wrap(~factor(gene), ncol=5, scales="free")+
   theme_bw()+
   theme(strip.text=element_text(size=12),
         axis.title.y=element_blank(),
         ##axis.text.x=element_text(size=8),
         axis.text.y=element_markdown())
         ##plot.margin=unit(c(0.1, 0.5, 0.1, 0.1), "cm"))
   
figfn <- paste(outdir, "Figure4.3_", ii, "_forestplot.png", sep="")
png(filename=figfn, width=1150, height=500, res=120)  
print(p3)
dev.off()



###############
### alcohol ###
###############

###
### category alcohol
ii <- "alcohol"
geneSel <- df%>%filter(SYMBOL%in%DEG_unq, category==ii)%>%pull(SYMBOL)%>%unique()

DEG_top <- res%>%dplyr::filter(gene%in%geneSel)%>%
    group_by(gene)%>%
    slice_max(order_by=abs(statistic), n=1)%>%ungroup()%>%
    arrange(desc(abs(statistic)))%>%pull(gene)%>%unique()



###
res2 <- res%>%filter(gene%in%DEG_top[1:5])%>%
   dplyr::select(gene, beta=estimate, SE=stderror, p.adjusted, celltype)

###
old2 <- old%>%filter(gene%in%DEG_top[1:5])%>%mutate(celltype="Bulk")%>%
   dplyr::select(gene, beta=log2FoldChange, SE=lfcSE, p.adjusted=padj, celltype)

 

###
plotDF <- rbind(res2, old2)
plotDF2 <- my_plotDF(plotDF)

plotDF2 <- plotDF2%>%
   mutate(p.adjusted=ifelse(is.na(p.adjusted), 1, p.adjusted),
  celltypeNew=gsub("_", " ", celltype),
  celltypeNew=ifelse(p.adjusted<0.1&abs(beta)>0, glue("<b style='font-size:10pt'><sup>*<b>{celltypeNew}"),
                     glue("</i>{celltypeNew}")),
          MCls_value=as.numeric(MCls_val[as.character(celltype)]),
          celltypeNew2=fct_reorder(celltypeNew, MCls_value),
          beta_upper=beta+1.96*SE,
          beta_lower=beta-1.96*SE)


###
### main
p3 <- ggplot(plotDF2, aes(x=beta, y=factor(celltypeNew2), color=factor(celltype)))+
   geom_point(shape=19, size=1.5, na.rm=T)+    
   geom_errorbarh(aes(xmax=beta_upper, xmin=beta_lower), size=0.5, height=0.2, na.rm=T)+
   scale_colour_manual(values=col_MCls, guide="none")+
   xlab("Log2FoldChange")+
   ##scale_y_discrete(breaks=c("ODC", "OPC", "Astrocyte", "Microglia", "DA", "Non_DA", "Pericyte", "Bulk"))+ 
   geom_vline(xintercept=0, size=0.25, linetype="dashed")+
   facet_wrap(~factor(gene), ncol=5, scales="free")+
   theme_bw()+
   theme(strip.text=element_text(size=12),
         axis.title.y=element_blank(),
         ##axis.text.x=element_text(size=8),
         axis.text.y=element_markdown())
         ##plot.margin=unit(c(0.1, 0.5, 0.1, 0.1), "cm"))
  
figfn <- paste(outdir, "Figure4.3_", ii, "_forestplot.png", sep="")
png(filename=figfn, width=950, height=250, res=120)  
print(p3)
dev.off()


### source plots data
x1 <- plotDF2%>%
   mutate(p.adjusted=ifelse(is.na(p.adjusted), 1, p.adjusted),
  celltypeNew=gsub("_", " ", celltype),
    MCls_value=as.numeric(MCls_val[as.character(celltype)]),
          beta_upper=beta+1.96*SE,
          beta_lower=beta-1.96*SE)

x1$traits <- "Alcohol"



#############
### smoke ### 
#############


###
### category smoke
ii <- "smoke"
geneSel <- df%>%filter(SYMBOL%in%DEG_unq, category==ii)%>%pull(SYMBOL)%>%unique()

DEG_top <- res%>%dplyr::filter(gene%in%geneSel)%>%
    group_by(gene)%>%
    slice_max(order_by=abs(statistic), n=1)%>%ungroup()%>%
    arrange(desc(abs(statistic)))%>%pull(gene)%>%unique()



###
res2 <- res%>%filter(gene%in%DEG_top[1:5])%>%
   dplyr::select(gene, beta=estimate, SE=stderror, p.adjusted, celltype)

###
old2 <- old%>%filter(gene%in%DEG_top[1:5])%>%mutate(celltype="Bulk")%>%
   dplyr::select(gene, beta=log2FoldChange, SE=lfcSE, p.adjusted=padj, celltype)

 

###
plotDF <- rbind(res2, old2)
plotDF2 <- my_plotDF(plotDF)

plotDF2 <- plotDF2%>%
   mutate(p.adjusted=ifelse(is.na(p.adjusted), 1, p.adjusted),
  celltypeNew=gsub("_", " ", celltype),
  celltypeNew=ifelse(p.adjusted<0.1&abs(beta)>0, glue("<b style='font-size:10pt'><sup>*<b>{celltypeNew}"),
                     glue("</i>{celltypeNew}")),
          MCls_value=as.numeric(MCls_val[as.character(celltype)]),
          celltypeNew2=fct_reorder(celltypeNew, MCls_value),
          beta_upper=beta+1.96*SE,
          beta_lower=beta-1.96*SE)


###
### main
p3 <- ggplot(plotDF2, aes(x=beta, y=factor(celltypeNew2), color=factor(celltype)))+
   geom_point(shape=19, size=1.5, na.rm=T)+    
   geom_errorbarh(aes(xmax=beta_upper, xmin=beta_lower), size=0.5, height=0.2, na.rm=T)+
   scale_colour_manual(values=col_MCls, guide="none")+
   xlab("Log2FoldChange")+
   ##scale_y_discrete(breaks=c("ODC", "OPC", "Astrocyte", "Microglia", "DA", "Non_DA", "Pericyte", "Bulk"))+ 
   geom_vline(xintercept=0, size=0.25, linetype="dashed")+
   facet_wrap(~factor(gene), ncol=5, scales="free")+
   theme_bw()+
   theme(strip.text=element_text(size=12),
         axis.title.y=element_blank(),
         ##axis.text.x=element_text(size=8),
         axis.text.y=element_markdown())
         ##plot.margin=unit(c(0.1, 0.5, 0.1, 0.1), "cm"))
  
figfn <- paste(outdir, "Figure4.3_", ii, "_forestplot.png", sep="")
png(filename=figfn, width=950, height=250, res=120)  
print(p3)
dev.off()


### source data

x2 <- plotDF2%>%
   mutate(p.adjusted=ifelse(is.na(p.adjusted), 1, p.adjusted),
  celltypeNew=gsub("_", " ", celltype),
 ## celltypeNew=ifelse(p.adjusted<0.1&abs(beta)>0, glue("<b style='font-size:10pt'><sup>*<b>{celltypeNew}"),
##                     glue("</i>{celltypeNew}")),
          MCls_value=as.numeric(MCls_val[as.character(celltype)]),
 ##         celltypeNew2=fct_reorder(celltypeNew, MCls_value),
          beta_upper=beta+1.96*SE,
          beta_lower=beta-1.96*SE)
x2$traits <- "Smoking"





################
### caffeine ### 
################


###
### category smoke
ii <- "caffeine"
geneSel <- df%>%filter(SYMBOL%in%DEG_unq, category==ii)%>%pull(SYMBOL)%>%unique()

DEG_top <- res%>%dplyr::filter(gene%in%geneSel)%>%
    group_by(gene)%>%
    slice_max(order_by=abs(statistic), n=1)%>%ungroup()%>%
    arrange(desc(abs(statistic)))%>%pull(gene)%>%unique()



###
res2 <- res%>%filter(gene%in%DEG_top[1:5])%>%
   dplyr::select(gene, beta=estimate, SE=stderror, p.adjusted, celltype)

###
old2 <- old%>%filter(gene%in%DEG_top[1:5])%>%mutate(celltype="Bulk")%>%
   dplyr::select(gene, beta=log2FoldChange, SE=lfcSE, p.adjusted=padj, celltype)

 

###
plotDF <- rbind(res2, old2)
plotDF2 <- my_plotDF(plotDF)

plotDF2 <- plotDF2%>%
   mutate(p.adjusted=ifelse(is.na(p.adjusted), 1, p.adjusted),
  celltypeNew=gsub("_", " ", celltype),
  celltypeNew=ifelse(p.adjusted<0.1&abs(beta)>0, glue("<b style='font-size:10pt'><sup>*<b>{celltypeNew}"),
                     glue("</i>{celltypeNew}")),
          MCls_value=as.numeric(MCls_val[as.character(celltype)]),
          celltypeNew2=fct_reorder(celltypeNew, MCls_value),
          beta_upper=beta+1.96*SE,
          beta_lower=beta-1.96*SE)


###
### main
p3 <- ggplot(plotDF2, aes(x=beta, y=factor(celltypeNew2), color=factor(celltype)))+
   geom_point(shape=19, size=1.5, na.rm=T)+    
   geom_errorbarh(aes(xmax=beta_upper, xmin=beta_lower), size=0.5, height=0.2, na.rm=T)+
   scale_colour_manual(values=col_MCls, guide="none")+
   xlab("Log2FoldChange")+
   ##scale_y_discrete(breaks=c("ODC", "OPC", "Astrocyte", "Microglia", "DA", "Non_DA", "Pericyte", "Bulk"))+ 
   geom_vline(xintercept=0, size=0.25, linetype="dashed")+
   facet_wrap(~factor(gene), ncol=5, scales="free")+
   theme_bw()+
   theme(strip.text=element_text(size=12),
         axis.title.y=element_blank(),
         ##axis.text.x=element_text(size=8),
         axis.text.y=element_markdown())
         ##plot.margin=unit(c(0.1, 0.5, 0.1, 0.1), "cm"))
  
figfn <- paste(outdir, "Figure4.3_", ii, "_forestplot2.png", sep="")
png(filename=figfn, width=950, height=250, res=120)  
print(p3)
dev.off()



x3 <- plotDF2%>%
   mutate(p.adjusted=ifelse(is.na(p.adjusted), 1, p.adjusted),
  celltypeNew=gsub("_", " ", celltype),
  ###celltypeNew=ifelse(p.adjusted<0.1&abs(beta)>0, glue("<b style='font-size:10pt'><sup>*<b>{celltypeNew}"),
  ##                   glue("</i>{celltypeNew}")),
          MCls_value=as.numeric(MCls_val[as.character(celltype)]),
  ##        celltypeNew2=fct_reorder(celltypeNew, MCls_value),
          beta_upper=beta+1.96*SE,
          beta_lower=beta-1.96*SE)

x3$traits <- "Caffeine"


xcomb <- rbind(x1, x2, x3)

opfn <- paste(outdir, "TableS_Figure5_forestplot.txt", sep="")
write.table(xcomb, opfn, row.names=F, col.names=T, quote=F, sep="\t")
