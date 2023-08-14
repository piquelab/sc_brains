##
library(Matrix)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(qvalue)

##
library(DESeq2)
library(annotables)
library(biobroom)
library(Seurat)
library(cowplot)
library(ggrepel)

library(openxlsx)
##
rm(list=ls())

outdir <- "./2_TWAS_DEG.outs/Response_to_reviews/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


###
### Twas Traits
### Last modified Jun-2-2023, by Julong wei


##############################
### read differential data ###
##############################



## old
## col_MCls <- c("ODC"="#fb8072", "Astrocyte"="#b3de69", "Microglia"="#bebada", "OPC"="#8dd3c7",
##    "Non_DA"="#fdb462",  "Pericyte"="#ffffb3", "DA"="#80b1d3", "Endothelial"="#fccde5",
##     "T-cell"="#35978f", "Ependymal"="#828282", "Union"="#ca0020")
###
MCls_name <- c("0"="ODC", "1"="Astrocyte", "2"="Microglia", "3"="OPC",
     "4"="Non_DA", "6"="Pericyte", "7"="DA", "8"="Endothelial", "9"="T-cell", "12"="Ependymal")
###
MCls_val <- c("ODC"=1, "OPC"=2, "Astrocyte"=3, "Microglia"=4, "DA"=5, "Non_DA"=6, 
    "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10, "Union"=11, "Old_Bulk"=12)
###
col_MCls <- c("ODC"="#fb8072", "OPC"="#35978f", "Astrocyte"="#92cd2d", "Microglia"="#c38dc4", 
   "DA"="#80b1d3", "Non_DA"="#fc9016",  "Pericyte"="#bebada",  "Endothelial"="#fccde5",
    "T-cell"="#8dd3c7",  "Ependymal"="#828282", "Union"="grey70", "Old_Bulk"="grey70")



### 
## Differential results
fn <- "../Merge-2022-06-02/4.0_Diff.cluster.outs/Correct_results/Filter2_default/1_DESeq_results_default.rds"
resDiff <- read_rds(fn)%>%drop_na(p.value)
resDiff2 <- resDiff%>%filter(p.adjusted<0.1, abs(estimate)>0.25)

celltypes <- sort(unique(resDiff2$celltype))

geneall <- unique(resDiff$gene)



## resMat <- resDiff%>%pivot_wider(id_cols=SYMBOL, names_from=MCls, values_from=is_DEG, values_fill=0)%>%as.data.frame()

## resMat2 <- resMat%>%mutate(Cluster7=rowSums(resMat[,-1]), Cluster7=ifelse(Cluster7>0, 1, 0))


##############################
### summary total overlap ####
##############################

###
### GTEx 
twas_new <- read.xlsx("TWAS_SUD_nmh_2023.xlsx")
twas_new <- twas_new%>%dplyr::select(gene, SYMBOL=gene_name, pval=pvalue)%>%
    mutate(gene=gsub("\\..*", "", gene), FDR=qvalue(pval)$qvalues)
twas_new2 <- twas_new%>%dplyr::filter(FDR<0.1)%>%mutate(traits="addiction-rf", category="addiction-GTEx")%>%
   dplyr::select(gene, pval, FDR, traits, category, SYMBOL)

###
opfn <- "TWAS_SUD_nmh_sig.xlsx"
write.xlsx(twas_new2, file=opfn, overwrite=T)


###
### PsychENCODE dataset
twas_new <- read.xlsx("TWAS_SUD_PsychENCODE_nmh_2023.xlsx")%>%drop_na(zscore, pvalue)

twas_new <- twas_new%>%dplyr::select(gene, zscore, pval=pvalue)%>%
    mutate(gene=gsub("\\..*", "", gene), FDR=qvalue(pval)$qvalues)
 
twas_new2 <- twas_new%>%
    dplyr::filter(FDR<0.1)%>%
    mutate(traits="addiction-rf", category="addiction-PsychoENCODE")

geneAnno <- bitr(twas_new2$gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)

twas_new2 <- twas_new2%>%left_join(geneAnno, by=c("gene"="ENSEMBL"))%>%drop_na(SYMBOL)
###
opfn <- "TWAS_SUD_nmh_PsychENCODE_sig.xlsx"
write.xlsx(twas_new2, opfn, overwrite=T)




           
###
### summary total overlap

gene_comb <- openxlsx::read.xlsx("./1_TWAS.outs/1_TWAS_traits.xlsx") 
trait_DF <- openxlsx::read.xlsx("./1_TWAS.outs/traits.xlsx")%>%dplyr::select(traits, Relevant)
gene_comb2 <- gene_comb%>%left_join(trait_DF, by="traits")%>%filter(Relevant==1)
gene_comb2 <- gene_comb2%>%dplyr::select(-Relevant)%>%filter(!category%in%c("ADHD", "parkinson"))


x3 <- gene_comb2%>%filter(SYMBOL%in%geneall)

###
x1 <- unique(gene_comb2$SYMBOL)
x2 <- unique(read.xlsx("TWAS_SUD_nmh_sig.xlsx")$SYMBOL)
x3 <- unique(read.xlsx("TWAS_SUD_nmh_PsychENCODE_sig.xlsx")$SYMBOL)

twas_gene <- intersect(c(x1, x2, x3), geneall)
DEG <- resDiff2%>%pull(gene)%>%unique()
olap <- intersect(DEG, twas_gene)




#####################################
### TWAS genes  for each category ###
#####################################



gene_comb <- openxlsx::read.xlsx("./1_TWAS.outs/1_TWAS_traits.xlsx") 
trait_DF <- openxlsx::read.xlsx("./1_TWAS.outs/traits.xlsx")%>%dplyr::select(traits, Relevant)
gene_comb2 <- gene_comb%>%left_join(trait_DF, by="traits")%>%filter(Relevant==1)
gene_comb2 <- gene_comb2%>%dplyr::select(-Relevant)%>%filter(!category%in%c("ADHD", "parkinson"))

gene_comb2 <- gene_comb2%>%dplyr::select(SYMBOL, category, traits)


### NEW SUD paper
x2 <- read.xlsx("TWAS_SUD_nmh_sig.xlsx")%>%dplyr::select(SYMBOL, category=traits)
x2$traits <- "addiction-rf"
x3 <- read.xlsx("TWAS_SUD_nmh_PsychENCODE_sig.xlsx")%>%dplyr::select(SYMBOL, category=traits, traits)
x3$traits <- "addiction-rf"

## contrast traits
fn <- "./1_TWAS.outs/1.2_TWAS_traits.xlsx"
x4 <- openxlsx::read.xlsx(fn)%>%dplyr::select(SYMBOL, category, traits)


##
gene_comb2 <- rbind(gene_comb2, x2, x3, x4)

gene_comb3 <- gene_comb2%>%filter(SYMBOL%in%geneall)


###
### trait of interest, supplemental tables

## gene_comb <- openxlsx::read.xlsx("./1_TWAS.outs/1_TWAS_traits.xlsx") 
## trait_DF <- openxlsx::read.xlsx("./1_TWAS.outs/traits.xlsx")%>%dplyr::select(traits, Relevant)
## gene_comb2 <- gene_comb%>%left_join(trait_DF, by="traits")%>%filter(Relevant==1)
## gene_comb2 <- gene_comb2%>%dplyr::select(-Relevant)%>%filter(!category%in%c("ADHD", "parkinson"))

## ##
## x2 <- read.xlsx("TWAS_SUD_nmh_sig.xlsx")%>%mutate(zscore=NA)%>%
##     dplyr::select(gene, zscore, pval, FDR, traits, category, SYMBOL)
## x2$category <- "addiction-GTEx"

## x3 <- read.xlsx("TWAS_SUD_nmh_PsychENCODE_sig.xlsx")

## xcomb <- rbind(gene_comb2, x2, x3)

## xcomb2 <- xcomb%>%mutate(category=gsub("-.*", "", category))
## summ <- xcomb2%>%group_by(category, traits)%>%summarise(nsig_twas=length(unique(SYMBOL)), .groups="drop")
## fn <- paste(outdir, "TableS_traits_list.xlsx", sep="")
## write.xlsx(summ, fn, overwrite=T)

## ###
## ### gene list for each trait
## fn <- paste(outdir, "TableS_genes_SUD.xlsx", sep="")
## write.xlsx(xcomb, fn, overwrite=T)




####################
### overlap DEGs ###
####################
gene_comb3 <- gene_comb3%>%mutate(category=gsub("-rf", "", category))

category <- sort(unique(gene_comb3$category))
DF <- NULL
for (oneMCl in celltypes){
    ##
    DEG <- resDiff2%>%filter(celltype==oneMCl)%>%pull(gene)%>%unique()
    ##
    df2 <- map_dfr(category, function(ii){
       ##
       twas_gene <- gene_comb3%>%filter(category==ii)%>%pull(SYMBOL)
       nolap <- length(intersect(DEG, twas_gene)) 
       df2 <- data.frame(MCls=oneMCl, category=ii, nolap=nolap, nDEG=length(DEG))
    })
    
    ##union_twas_traits
    twas_gene <- unique(gene_comb3$SYMBOL)
    nolap <- length(intersect(DEG, twas_gene))    
    df3 <- data.frame(MCls=oneMCl, category="Union_twas", nolap=nolap, nDEG=length(DEG))
    DF <- rbind(DF, df2, df3)
}


###
### union of DEGs across cell-types
DEG <- resDiff2%>%pull(gene)%>%unique()
##
df2 <- map_dfr(category, function(ii){
##
   twas_gene <- gene_comb3%>%filter(category==ii)%>%pull(SYMBOL)
   nolap <- length(intersect(DEG, twas_gene)) 
   df2 <- data.frame(MCls="Union", category=ii, nolap=nolap, nDEG=length(DEG))
})

##union_twas_traits
twas_gene <- unique(gene_comb3$SYMBOL)
nolap <- length(intersect(DEG, twas_gene))    
df3 <- data.frame(MCls="Union", category="Union_twas", nolap=nolap, nDEG=length(DEG))
##
DF <- rbind(DF, df2, df3)



###
### old bulk data
fn <- "/wsu/home/groups/bannonlab/manal/DESeq2_WGCNA_DownstreamAnalysis/Final_Runs/DESEQ2_TEST_3-28-18/DESEQ2_IHW-5-28-18.csv"
old <- read.csv(fn)%>%drop_na(pvalue)%>%dplyr::rename("gene"="Associated.Gene.Name")
old_DEG <- old%>%filter(padj<0.1)%>%pull(gene)%>%unique()

df2 <- map_dfr(category, function(ii){
##
   twas_gene <- gene_comb3%>%filter(category==ii)%>%pull(SYMBOL)
   nolap <- length(intersect(old_DEG, twas_gene)) 
   df2 <- data.frame(MCls="Old_Bulk", category=ii, nolap=nolap, nDEG=length(old_DEG))
})

##union_twas_traits
twas_gene <- unique(gene_comb3$SYMBOL)
nolap <- length(intersect(old_DEG, twas_gene))    
df3 <- data.frame(MCls="Old_Bulk", category="Union_twas", nolap=nolap, nDEG=length(DEG))
DF <- rbind(DF, df2, df3)




### final plot data
plotDF <- DF%>%
    mutate(percent=nolap/nDEG, percent2=ifelse(nolap==0, NA, percent),
    MCl_value=MCls_val[as.character(MCls)],
    MCl2=fct_reorder(MCls, MCl_value))

 
###
plotDF2 <- plotDF%>%filter(!category%in%c("Union_twas", "height", "CAD", "BMI"))
p <- ggplot(plotDF2, aes(x=MCl2, y=category, size=percent2, color=MCls))+
   geom_point(alpha=0.8)+
   scale_color_manual(values=col_MCls, guide="none")+ 
   ## scale_color_manual(values=c("addiction"="#e41a1c", "ADHD"="#377eb8", "alcohol"="#4daf4a",
   ##    "caffeine"="#984ea3", "parkinson"="#ff7f00", "smoke"="#e41a1c"),
    ##    guide="none")+
   scale_size_binned("Percent", breaks=waiver(), n.breaks=4,
      guide=guide_bins(show.limits=T, axis=T,
          override.aes=list(color="#fb8072"),             
          axis.show=arrow(length=unit(1.5, "mm"), ends="both"),
          keywidth=grid::unit(0.2, "lines"),
          keyheight=grid::unit(0.6, "lines") ))+ 
   scale_x_discrete(labels=c("ODC"="ODC", "OPC"="OPC", "Astrocyte"="Astrocyte",
       "Microglia"="Microglia", "Non_DA"="Non DA", "Pericyte"="Pericyte", "Union"="Union", "Old_Bulk"="Bulk"))+
   theme_bw()+
   theme(axis.title=element_blank(),
         axis.text.x=element_text(angle=60, hjust=1, size=10),
         axis.text.y=element_text(size=10))
 
### 
figfn <- paste(outdir, "Figure1_bubble.png", sep="")
png(figfn, width=620, height=400, res=120)
print(p)
dev.off()






#################################################
### Table of overlapping DEGs with TWAS genes ###
#################################################

## summary #traits and #genes
summ_twas <- gene_comb3%>%group_by(category)%>%
    summarise(ntraits=length(unique(traits)), ngene=length(unique(SYMBOL)), .groups="drop")%>%ungroup()%>%
    arrange(category)
tmp <- data.frame(category="Union_twas", ntraits=length(unique(gene_comb3$traits)), ngene=length(unique(gene_comb3$SYMBOL)))
###
summ_twas <- rbind(summ_twas, tmp)



summ <- DF%>%arrange(category)%>%pivot_wider(id_cols=category, names_from=MCls, values_from=nolap)
summ2 <- summ[,c("category", "ODC", "OPC", "Astrocyte", "Microglia", "Non_DA", "Pericyte", "Union", "Old_Bulk")]

comb <- cbind(summ_twas, summ2[,-1])
##
opfn <- paste(outdir, "1_summary_olap.xlsx", sep="")
openxlsx::write.xlsx(comb, opfn, overwrite=T)







#######################
### proportion test ###
#######################


celltypes <- c("ODC", "OPC", "Astrocyte", "Microglia", "Non_DA", "Pericyte", "Union", "Old_Bulk") 
traits <- c("addiction-rf", "alcohol", "BMI", "CAD", "caffeine", "height",  "marijuana", "smoke")


###
### proportion test, which trait are enriched for each cell-type
res_prop <- map_dfc(celltypes, function(oneMCl){
   ###
   pval <- sapply(traits, function(ii){
      ##
      nolap_obs <- DF%>%filter(MCls==oneMCl, category==ii)%>%pull(nolap)
      nolap_exp <- DF%>%filter(MCls==oneMCl, category=="Union_twas")%>%pull(nolap)
      ntwas_each <- summ_twas%>%filter(category==ii)%>%pull(ngene)
      ntwas_total <- summ_twas%>%filter(category=="Union_twas")%>%pull(ngene)
      ###
      res <- prop.test(c(nolap_obs, nolap_exp), c(ntwas_each, ntwas_total), alternative="greater")
      pval <- res$p.value       
   })
   pval
})

res_prop <- as.data.frame(res_prop)
names(res_prop) <- celltypes
rownames(res_prop) <- traits




#####################################
### if enriched in some cell-types ##
#####################################

DEG_union <- resDiff2%>%pull(gene)%>%unique()
DF_deg <- resDiff2%>%group_by(celltype)%>%summarise(ngene=length(unique(gene)), .groups="drop")%>%ungroup()

###
### proportion test, which trait are enriched for each cell-type
res_prop <- map_dfc(celltypes[1:6], function(oneMCl){
   ###
   pval <- sapply(traits, function(ii){
      ##
      nolap_obs <- DF%>%filter(MCls==oneMCl, category==ii)%>%pull(nolap)
      nolap_exp <- DF%>%filter(MCls=="Union", category==ii)%>%pull(nolap)
      ndeg_each <- DF_deg%>%filter(celltype==oneMCl)%>%pull(ngene)
      ndeg_total <- length(DEG_union)
      ###
      res <- prop.test(c(nolap_obs, nolap_exp), c(ndeg_each, ndeg_total), alternative="greater")
      pval <- res$p.value       
   })
   pval
})

res_prop <- as.data.frame(res_prop)
names(res_prop) <- celltypes[1:6]
rownames(res_prop) <- traits






#####################################################################
### Barplot of number of TWAS-gene overlap DEG for each cell type ###
#####################################################################


gene_comb <- openxlsx::read.xlsx("./1_TWAS.outs/1_TWAS_traits.xlsx") 
trait_DF <- openxlsx::read.xlsx("./1_TWAS.outs/traits.xlsx")%>%dplyr::select(traits, Relevant)
gene_comb2 <- gene_comb%>%left_join(trait_DF, by="traits")%>%filter(Relevant==1)
gene_comb2 <- gene_comb2%>%dplyr::select(-Relevant)%>%filter(!category%in%c("ADHD", "parkinson"))
x1 <- unique(gene_comb2$SYMBOL)
x2 <- unique(read.xlsx("TWAS_SUD_nmh_sig.xlsx")$SYMBOL)
x3 <- unique(read.xlsx("TWAS_SUD_nmh_PsychENCODE_sig.xlsx")$SYMBOL)

twas_gene <- intersect(c(x1, x2, x3), geneall)

### 
###
celltype2 <- sort(unique(resDiff2$celltype))
plotDF <- NULL
for (ii in celltype2){
###
   DEG <- resDiff2%>%filter(celltype==ii)%>%pull(gene)%>%unique()
   x <- data.frame(celltype=ii, nolap=length(intersect(DEG, twas_gene)))
   plotDF <- rbind(plotDF, x)
}

plotDF2 <- plotDF%>%
    mutate(MCl_value=as.numeric(MCls_val[as.character(celltype)]),
    MCl2=fct_reorder(celltype, MCl_value))

###
### barplots 
###plotDF2 <- plotDF2%>%filter(!celltype%in%c("Union", "Old_Bulk"))
p <- ggplot(plotDF2, aes(x=nolap, y=MCl2, fill=MCl2))+
    geom_bar(stat="identity")+
    scale_fill_manual(values=col_MCls, guide="none")+
    ggtitle("#DEGs overlapping TWAS SUD")+
    scale_y_discrete(labels=c("ODC"="ODC", "OPC"="OPC", "Astrocyte"="Astrocyte",
       "Microglia"="Microglia", "Non_DA"="Non DA", "Pericyte"="Pericyte"))+    
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(size=12),
          plot.title=element_text(hjust=0.5, size=12))

figfn <- paste(outdir, "Figure2_olap_barplot.png", sep="")
png(figfn, width=450, height=500, res=120)
print(p)
dev.off()

 

DEGs <- resDiff2%>%pull(gene)%>%unique()
nolap <- length(intersect(DEGs, twas_gene))

 


####################################
### scatter plots of DEG vs TWAS ###
####################################

### 
## Differential results
fn <- "../Merge-2022-06-02/4.0_Diff.cluster.outs/Correct_results/Filter2_default/1_DESeq_results_default.rds"
resDiff <- read_rds(fn)%>%drop_na(p.value)
resDiff2 <- resDiff%>%filter(p.adjusted<0.1, abs(estimate)>0.25)

celltypes <- sort(unique(resDiff2$celltype))

geneall <- unique(resDiff$gene)


###
### gwas results

twas <- read.xlsx("TWAS_SUD_nmh_PsychENCODE_sig.xlsx")
twas2 <- twas%>%dplyr::select(gene=SYMBOL, zscore_TWAS=zscore)


celltypes <- c("ODC", "OPC", "Astrocyte", "Microglia")
plotDF <- map_dfr(celltypes, function(ii){
    ##
    x <- resDiff2%>%filter(celltype==ii)
    x <- x%>%dplyr::select(gene, zscore_DE=statistic, celltype)
    x <- x%>%inner_join(twas2, by="gene")
    cat(ii, nrow(x), "\n")
    x
})    

##
plotDF2 <- plotDF%>%
    mutate(MCls_value=as.numeric(MCls_val[celltype]),
        celltype2=fct_reorder(celltype, MCls_value))

p2 <- ggplot(plotDF2, aes(x=zscore_DE, y=zscore_TWAS, color=factor(celltype2)))+
   geom_point(shape=24)+
   geom_text_repel(aes(x=zscore_DE, y=zscore_TWAS, label=gene),
                   size=2, fontface="italic", max.overlaps=20)+
   scale_color_manual(values=col_MCls, guide="none")+
   xlab("zscore of DE")+
   ylab("zscore of addiction-rf GWAS")+
   geom_hline(yintercept=0, linetype="dashed", color="grey30")+
   geom_vline(xintercept=0, linetype="dashed", color="grey30")+
   facet_wrap(~celltype2, ncol=2, scales="fixed")+
   theme_bw()+
   theme(strip.text=element_text(size=12),
         axis.title=element_text(size=10))

figfn <- paste(outdir, "Figure3_TWAS_scatter.png", sep="")
png(figfn, width=620, height=600, res=120)
print(p2)
dev.off()



### End










## ####
## ####
## ### current color
## col_MCls <- c("ODC"="#fb8072", "Astrocyte"="#92cd2d", "Microglia"="#c38dc4", "OPC"="#35978f",
##    "Non_DA"="#fc9016",  "Pericyte"="#bebada", "DA"="#80b1d3", "Endothelial"="#fccde5",
##     "T-cell"="#8dd3c7",  "Ependymal"="#828282")

## ###
## MCls_name <- c("0"="ODC", "1"="Astrocyte", "2"="Microglia", "3"="OPC",
##      "4"="Non_DA", "6"="Pericyte", "7"="DA", "8"="Endothelial", "9"="T-cell", "12"="Ependymal")
## ###
## MCls_val <- c("ODC"=1, "OPC"=2, "Astrocyte"=3, "Microglia"=4, "DA"=5, "Non_DA"=6, 
##     "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)

## ###
## ### Differential results
## fn <- "../Merge-2022-06-02/4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds"
## resDiff <- read_rds(fn)%>%drop_na(p.value)%>%
##     mutate(celltype=MCls_name[as.character(MCls)])%>%
##     filter(p.adjusted<0.1, abs(estimate)>0.25)


## ###
## gene_comb <- openxlsx::read.xlsx("./1_TWAS.outs/1_TWAS_traits.xlsx") 
## trait_DF <- openxlsx::read.xlsx("./1_TWAS.outs/traits.xlsx")%>%dplyr::select(traits, Relevant)
## gene_comb2 <- gene_comb%>%left_join(trait_DF, by="traits")%>%filter(Relevant==1)



## ###
## ###
## fn <- "./spredixcan/spredixcan-Brain_Substantia_nigra.tsv.gz"
## twas <- fread(fn, data.table=F)
## twas <- twas%>%distinct(gene_name, .keep_all=T)%>%column_to_rownames(var="gene_name")

## ###
## traits <- colnames(twas)
## gene <- rownames(twas)
## geneAnno <- bitr(gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)






## ######################
## #### scatter plots ###
## ######################

## fn <- "../Merge-2022-06-02/4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds"
## resDiff <- read_rds(fn)%>%drop_na(p.value)%>%
##     mutate(celltype=MCls_name[as.character(MCls)])

## celltypes <- sort(unique(resDiff$celltype))

## ## resMat <- resDiff%>%pivot_wider(id_cols=SYMBOL, names_from=MCls, values_from=is_DEG, values_fill=0)%>%as.data.frame()

## ## resMat2 <- resMat%>%mutate(Cluster7=rowSums(resMat[,-1]), Cluster7=ifelse(Cluster7>0, 1, 0))



## trait_DF <- openxlsx::read.xlsx("./1_TWAS.outs/traits.xlsx")%>%dplyr::select(traits, Relevant)%>%filter(Relevant==1)

## ###
## ###

## fn <- "./spredixcan/spredixcan-Brain_Substantia_nigra.tsv.gz"
## twas <- fread(fn, data.table=F)
## twas <- twas%>%distinct(gene_name, .keep_all=T) ##%>%column_to_rownames(var="gene_name")

## ###
## traits <- colnames(twas)
## gene <- twas$gene_name
## geneAnno <- bitr(gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)

## twas2 <- twas%>%inner_join(geneAnno, by=c("gene_name"="ENSEMBL"))

 
## ###
## plotDF <- twas2%>%dplyr::select(gene=SYMBOL, zscore_twas=PGC_ADHD_EUR_2017)

## resDiff2 <- resDiff%>%dplyr::select(gene, zscore_DE=statistic, celltype)

## plotDF <- plotDF%>%left_join(resDiff2, by="gene")


## plotDF2 <- plotDF%>%
##     mutate(MCls_value=as.numeric(MCls_val[as.character(celltype)]),
##         celltype2=fct_reorder(celltype, MCls_value))%>%drop_na(celltype)                   

    
## ###
## ####
## p <- ggplot(plotDF2, aes(x=zscore_twas, y=zscore_DE))+
##     geom_point(size=0.2, color="grey40")+
##     facet_wrap(~celltype2, ncol=4, scales="fixed")+
##     theme_bw()

## figfn <- paste(outdir, "Figure2_zscore_compare.png", sep="")
## png(figfn, width=880, height=500, res=120)
## print(p)
## dev.off()

    


