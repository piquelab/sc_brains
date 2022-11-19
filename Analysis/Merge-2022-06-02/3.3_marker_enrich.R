##
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)

## library(rmarkdown)
## library(knitr)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(ggtext, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(glue)

rm(list=ls())

outdir <- "./3_cluster.outs/plots_pub"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###
###

res <- read_rds("./3_cluster.outs/2.0_cluster14_allmarker.rds")
res2 <- res%>%filter(p_val_adj<0.1, avg_log2FC>0.25)

## res2%>%group_by(cluster)%>%summarise(ngene=length(unique(gene)), .groups="drop")

combs <- sort(unique(res$cluster))
results <- lapply(combs, function(ii){
   ### 
   cat(ii, "\n") 
   ## 
   BG <- res%>%filter(cluster==ii)%>%pull(gene)%>%unique()
   BG_df <- bitr(BG, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=org.Hs.eg.db)
   ##
    
   DEG <- res2%>%filter(cluster==ii)%>%pull(gene)%>%unique()
   DEG_df <- bitr(DEG, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=org.Hs.eg.db)  
 
   ego2 <- enrichGO(gene=DEG_df$ENTREZID,
                 universe=BG_df$ENTREZID,
                 OrgDb=org.Hs.eg.db,
                 ont="BP",
                 minGSSize=1,
                 maxGSSize=10000,
                 pvalueCutoff=1, qvalueCutoff=1, readable=T)
   ego2
})

names(results) <- combs

cg <- merge_result(results)

opfn <- "./3_cluster.outs/2_enriched.rds"
write_rds(cg, opfn)




#############
### plots ###
#############


cg <- read_rds("./3_cluster.outs/2_enriched.rds")
cg2 <- cg%>%mutate(maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))%>%
    filter(maxGSSize>=5&maxGSSize<=500, p.adjust<0.1)

###
MCls_name <- c("0"="0:ODC", "5"="5:ODC", "10"="10:ODC", "3"="3:OPC", "1"="1:Astrocyte",
   "2"="2:Microglia", "11"="11:Microglia", "7"="7:DA", "4"="4:Non_DA",
   "6"="6:Pericyte", "8"="8:Endothelial", "9"="9:T-cell", "12"="12:Ependymal")
###
MCls_value <- c("0"=0, "5"=1, "10"=2, "3"=3, "1"=4, "2"=5, "11"=6, "7"=7, "4"=8,
    "6"=9,  "8"=10, "9"=11, "12"=12)


###
cg2 <- cg2%>%mutate(ClusterNew=as.character(MCls_name[as.character(Cluster)]),
    CL_value=as.numeric(MCls_value[as.character(Cluster)]),
    ClusterNew2=fct_reorder(ClusterNew, CL_value))


p <- enrichplot::dotplot(cg2, x="ClusterNew2", showCategory=10)+
    scale_y_discrete(labels=function(y) str_wrap(y, width=60))+
    theme(axis.title=element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.text.y=element_text(size=12))
 
###
figfn <- paste(outdir, "/FigureS2_enrichGO.png", sep="")
png(figfn, width=1500, height=1200, res=120)
p
dev.off()
          
