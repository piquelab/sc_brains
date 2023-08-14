##
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)
library(qvalue)

## library(rmarkdown)
## library(knitr)
library(RColorBrewer)
library(cowplot)
library(ggplot2)
library(ggtext, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(glue)
library(openxlsx)

rm(list=ls())

outdir <- "./5_enriched.outs/Correct_results/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)




###
### differential results
## resDiff <- read_rds("./4.0_Diff.cluster.outs/Filter2_cells30_M2_batch/1.4_DESeq.results_CPM0.5_filter.rds")
## resDiff <- resDiff%>%drop_na(p.value)%>%
##     mutate(celltype=MCls_name[as.character(MCls)],
##        direction=ifelse(estimate>0, 1, 2), comb=paste(celltype, direction, sep="_"))%>%
##     as.data.frame()


fn <- "./4.0_Diff.cluster.outs/Correct_results/Filter2_default/1_DESeq_results_default.rds"
resDiff <- read_rds(fn)%>%drop_na(p.value)%>%
   mutate(direction=ifelse(estimate>0, 1, 2), comb=paste(celltype, direction, sep="_"))%>%
    as.data.frame()
resDiff2 <- resDiff%>%filter(p.adjusted<0.1, abs(estimate)>0.25)


###############################
### main, enrichGO analysis ###
###############################

combs <- sort(unique(resDiff2$comb))
results <- lapply(combs, function(ii){
   ###
   cat(ii, "\n") 
   oneMCl <- gsub("_[12]", "", ii)
   ## 
   BG <- resDiff%>%filter(celltype==oneMCl)%>%pull(gene)%>%unique()
   BG_df <- bitr(BG, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=org.Hs.eg.db)
   ##
    
   DEG <- resDiff2%>%filter(comb==ii)%>%pull(gene)%>%unique()
   DEG_df <- bitr(DEG, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=org.Hs.eg.db)  
 
   ego2 <- enrichGO(gene=DEG_df$ENTREZID,
                 universe=BG_df$ENTREZID,
                 OrgDb=org.Hs.eg.db,
                 ont="ALL",
                 minGSSize=1,
                 maxGSSize=10000,
                 pvalueCutoff=1, qvalueCutoff=1, readable=T)
   ego2
})

names(results) <- combs

cg <- merge_result(results)

opfn <- paste(outdir, "1_enriched.direction.rds", sep="")
write_rds(cg, opfn)


##########################################
### calcualte pvalue using fisher.test ###
##########################################
 
fn <- paste(outdir, "1_enriched.direction.rds", sep="")
res <- read_rds(fn)%>%
   mutate(DEG_in=as.numeric(gsub("/.*", "", GeneRatio)), DEG=as.numeric(gsub(".*/", "", GeneRatio)),
          BG_in=as.numeric(gsub("/.*", "", BgRatio)), BG=as.numeric(gsub(".*/", "", BgRatio)),
          DEG_notIN=DEG-DEG_in, notDEG_in=BG_in-DEG_in, notDEG_notIN=BG-DEG-notDEG_in)
x <- res@compareClusterResult
x2 <- x%>%dplyr::select(Cluster, ID, DEG_in, notDEG_in, DEG_notIN, notDEG_notIN)

### data frame
myFisher <- lapply(1:nrow(x2), function(i){
    ##
    dat <- matrix(as.numeric(x2[i, 3:6]), 2, 2)
    aux <- fisher.test(dat,alternative="greater")
    aux2 <- fisher.test(dat)
    ##
    aux2 <- data.frame(pval_mine=aux$p.value, odds_mine=aux2$estimate,
                       CI_lower=aux2$conf.int[1], CI_upper=aux2$conf.int[2])
    aux2
})
myFisher <- do.call(rbind, myFisher)
rownames(myFisher) <- NULL

x <- cbind(x, myFisher)

x <- x%>%group_by(ONTOLOGY)%>%mutate(p.adjust_across=p.adjust(pval_mine, method="BH"))%>%ungroup()


###
res@compareClusterResult <- x
opfn <- paste(outdir, "2_enriched_Fisher.test.rds", sep="")
write_rds(res, file=opfn)

##
### adjust p-value
## fn <- paste(outdir, "2_enriched_Fisher.test.rds", sep="")
## res <- read_rds(fn)
## ##
## x <- res@compareClusterResult
## x <- x%>%mutate(MCls=gsub("_[12]", "", Cluster))
## x <- x%>%group_by(MCls)%>%mutate(p.adjust_MCls=p.adjust(pval_mine, method="BH"))%>%ungroup()
 
## res@compareClusterResult <- x

## opfn <- paste(outdir, "2_enriched_Fisher.test.rds", sep="")
## write_rds(res, file=opfn)




###
### compared p-value

fn <- paste(outdir, "2_enriched_Fisher.test.rds", sep="")
res <- read_rds(fn)
x <- res@compareClusterResult

plotDF <- x%>%dplyr::select(Cluster, ID, x=pvalue, y=pval_mine)%>%
   mutate(MCls=gsub("_[12]", "", Cluster),
          dir2=ifelse(grepl("_1", Cluster), "Up", "Down"),
          Cluster2=paste(MCls, dir2, sep="_"),
          x=-log10(x), y=-log10(y))

p <- ggplot(plotDF, aes(x=x, y=y))+
   geom_point(size=0.5)+
   geom_abline(color="red")+
   xlab(bquote(~log[10]~"("~italic(p)~")"~"from ClusterProfiler"))+
   ylab(bquote(~log[10]~"("~italic(p)~")"~"from Fisher.test"))+ 
   facet_wrap(~Cluster2, nrow=3, ncol=4, scales="free")+
   theme_bw()+
   theme(axis.title=element_text(size=14),
         axis.text=element_text(size=12),
         strip.text=element_text(size=14))
 
figfn <- paste(outdir, "test/Figure0_compare_pval.png", sep="")
png(figfn, width=1000, height=800, res=120)
print(p)
dev.off()


###
###
## fn <- paste(outdir, "2_enriched_Fisher.test.rds", sep="")
## res <- read_rds(fn)
## x <- res@compareClusterResult
## x2 <- x%>%filter(p.adjust_mine<0.1, BG_in>=5&BG_in<=500, Cluster=="Non_DA_2")

## ###
## genes <- lapply(1:nrow(x2), function(i){
##    ##
##    gene <- str_split(x2$geneID[i], "/")
##    gene <- unlist(gene)
##    gene
## })
## genes <- unique(unlist(genes))






#############
### plots ###
#############

## MCls_val <- 1:12
## MCls <- c("ODC", "OPC", "Astrocyte", "Microglia",  "Non_DA", "Pericyte")
## names(MCls_val) <- c(paste(MCls, "_Down", sep=""), paste(MCls, "_Up", sep=""))

 
## cg <- read_rds("./5_enriched.outs/Correct_results/2_enriched_Fisher.test.rds")
## cg2 <-cg%>%dplyr::rename(p.adjust_across=p.adjust_mine)%>%
##     filter(BG_in>=5&BG_in<=500, p.adjust_across<0.2) ##, p.adjust<0.1|p.adjust_mine<0.1|p.adjust_MCls<0.1)

## ###
## cg2 <- cg2%>%mutate(MCls=gsub("_[12]", "", Cluster),
##     direction=gsub(".*_", "", Cluster),
##     dir2=ifelse(direction=="1", "Up", "Down"),
##     clusterNew=paste(MCls, dir2, sep="_"),
##     cluster_value=as.numeric(MCls_val[as.character(clusterNew)]),
##     clusterNew2=fct_reorder(clusterNew, cluster_value))

## p <- enrichplot::dotplot(cg2, x="clusterNew2", showCategory=5, color="p.adjust_across")+
##     scale_y_discrete(labels=function(y) str_wrap(y, width=80))+
##     theme(axis.title=element_blank(),
##           axis.text.x=element_text(angle=45, hjust=1, size=12),
##           axis.text.y=element_text(size=9))
 
## ###
## figfn <- paste(outdir, "Figure1_enrichGO.png", sep="")
## png(figfn, width=1500, height=1200, res=120)
## p
## dev.off()



###############
### summary ###
###############

###
### differential results
fn <- "./4.0_Diff.cluster.outs/Correct_results/Filter2_default/1_DESeq_results_default.rds"
resDiff <- read_rds(fn)%>%drop_na(p.value)%>%
   mutate(direction=ifelse(estimate>0, 1, 2), comb=paste(celltype, direction, sep="_"))%>%
    as.data.frame()
resDiff2 <- resDiff%>%filter(p.adjusted<0.1, abs(estimate)>0.25)
summDE <- resDiff2%>%group_by(comb)%>%summarise(ngene=n(),.groups="drop")
names(summDE) <- c("Cluster", "ngene")

###
### GO results
cg <- read_rds("./5_enriched.outs/Correct_results/2_enriched_Fisher.test.rds")
x <- cg@compareClusterResult
x2 <- x%>%filter(BG_in>=5, BG_in<=500, ONTOLOGY=="BP")

summ2 <- x2%>%group_by(Cluster)%>%
       summarise(ngo=n(), nsig_padj_across=sum(p.adjust_across<0.1))%>%ungroup()
##
summ2 <- summ2%>%
   mutate(MCls=gsub("_[12]", "", Cluster),
          Drt2=ifelse(grepl("_1", Cluster), "Up", "Down"),
          Cluster2=paste(MCls, Drt2, sep="_"))


###
### genes mapping to GO BP
comb <- sort(unique(x2$Cluster))
summBP <- NULL
for (ii in comb){
   ###
   cat(ii, "\n") 
   x3 <- x2%>%filter(Cluster==ii) 
   genes <- lapply(1:nrow(x3), function(i){
      ##
      gene <- str_split(x3$geneID[i], "/")
      gene <- unlist(gene)
      gene
   })
   genes <- unique(unlist(genes))
   summBP <- rbind(summBP, data.frame(Cluster=ii, ngene2=length(genes)))
}


### summary table
summ2 <- summ2%>%left_join(summDE, by="Cluster")%>%left_join(summBP, by="Cluster")

opfn <- paste(outdir, "0_summary.xlsx", sep="")
openxlsx::write.xlsx(summ2, file=opfn, overwrite=T)


############################################
### outputs of significantly enriched GO ###
############################################

cg <- read_rds("./5_enriched.outs/Correct_results/2_enriched_Fisher.test.rds")
x <- cg@compareClusterResult
x2 <- x%>%filter(BG_in>=5, BG_in<=500, ONTOLOGY=="BP", p.adjust_across<0.1)
###
###
x3 <- x2[,c(1:6, 10, 12:15, 19:23)]
x3 <- x3%>%
   mutate(MCls=gsub("_[12]", "", Cluster),
          direction=ifelse(grepl("_1", Cluster), "Up", "Down"),
          Cluster2=paste(MCls, direction, sep="_"))%>%
   dplyr::rename(pval=pval_mine, odds=odds_mine)
###
###


outdir2 <- "./5_enriched.outs/Correct_results/enriched_GO/"
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)

### output together

comb <- unique(x3$Cluster2)

xsort <- map_dfr(comb, ~x3%>%filter(Cluster2==.x)%>%arrange(pval))


opfn <- paste(outdir2, "combine_enriched_GO.xlsx", sep="")
openxlsx::write.xlsx(xsort, file=opfn, overwrite=T)


###
###

outdir2 <- "./5_enriched.outs/Correct_results/enriched_GO/"
fn <- paste(outdir2, "combine_enriched_GO.xlsx", sep="")
x <- read.xlsx(fn)

fn <- "./5_enriched.outs/Correct_results/1_enriched.direction.rds"
res <- read_rds(fn)








##################
### tree plots ###
##################

rm(list=ls())

outdir2 <- "./5_enriched.outs/Correct_results/treeplot/"
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)

fn <- "./4.0_Diff.cluster.outs/Correct_results/Filter2_default/1_DESeq_results_default.rds"
resDiff <- read_rds(fn)%>%drop_na(p.value)%>%
   mutate(direction=ifelse(estimate>0, 1, 2), comb=paste(celltype, direction, sep="_"))%>%
    as.data.frame()
resDiff2 <- resDiff%>%filter(p.adjusted<0.1, abs(estimate)>0.25)


###
###

cg <- read_rds("./5_enriched.outs/Correct_results/2_enriched_Fisher.test.rds")
x <- cg@compareClusterResult%>%filter(ONTOLOGY=="BP")


combs <- sort(unique(resDiff2$comb))
for ( ii in combs){
   ## 
   oneMCl <- gsub("_[12]", "", ii)
   dir2 <- ifelse(grepl("_1", ii), "Up", "Down")
    
   ## 
   BG <- resDiff%>%filter(celltype==oneMCl)%>%pull(gene)%>%unique()
   BG_df <- bitr(BG, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=org.Hs.eg.db)
   ##    
   DEG <- resDiff2%>%filter(comb==ii)%>%pull(gene)%>%unique()
   DEG_df <- bitr(DEG, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=org.Hs.eg.db)  

   ### enrichment procedure 
   cg2 <- enrichGO(gene=DEG_df$ENTREZID,
                 universe=BG_df$ENTREZID,
                 OrgDb=org.Hs.eg.db,
                 ont="BP",
                 minGSSize=1,
                 maxGSSize=10000,
                 pvalueCutoff=1, qvalueCutoff=1, readable=T) 

   ### 
   x2 <- x%>%filter(Cluster==ii)%>%dplyr::select(ID, p.adjust_across) 
   padj <- x2$p.adjust_across
   names(padj) <- as.character(x2$ID)
    

   ## ##
   ## res <- ego2@result
   ## res <- res%>%dplyr::left_join(x2, by="ID")
   ## ego2@result <- res
   cg3 <- cg2%>%mutate(GO_size=as.numeric(gsub("/.*", "", BgRatio)),
                       p.adjust=padj[as.character(ID)])%>%
       filter(GO_size>=5&GO_size<=500, padj<0.1)  ##%>%dplyr::select(-GO_size) , -p.adjust_across)

   ###
   ngo <- nrow(cg3)    
   if ( ngo>=10){
      ###
      cg3x <- pairwise_termsim(cg3) 
      p0 <- enrichplot::treeplot(cg3x, color="p.adjust", showCategory=50)+
          ggtitle(paste(oneMCl, dir2, sep="_"))+
          theme(plot.title=element_text(hjust=0.5, size=14))

      ### 
      figfn <- paste(outdir2, "Figure_", oneMCl, "_", dir2, "_treeplot.pdf", sep="") 
      if ( ngo>=30){ 
         pdf(figfn, width=10, height=10)
         print(p0)
         dev.off()
      }else{
         pdf(figfn, width=10, height=8)
         print(p0)
         dev.off()
      }   
   }else if (ngo>0){
      ## dot plots    
      p0 <- enrichplot::dotplot(cg3, color="p.adjust", showCategory=ngo)+
          ggtitle(paste(oneMCl, dir2, sep="_"))+
          theme(plot.title=element_text(hjust=0.5, size=14))
      ##  
      figfn <- paste(outdir2, "Figure_", oneMCl, "_", dir2, "_dotplot.pdf", sep="")
      pdf(figfn, width=8, height=8)
      print(p0)
      dev.off()
   }
    
   cat(ii, ngo, "\n") 
###
}





###
### outputs 

fn <- "./5_enriched.outs/Correct_results/2_enriched_Fisher.test.rds"
cg <- read_rds(fn)
x <- cg@compareClusterResult

###
### 
comb <- sort(unique(x$Cluster))
for (ii in comb){
   ##
   cat(ii, "\n") 
   direction <- ifelse(gsub(".*_", "", ii)=="1", "Up", "Down")
   MCl <- gsub("_[12]", "", ii)

   ## test GOs
   x2 <- x%>%filter(Cluster==ii)%>%arrange(pvalue)%>%
       dplyr::select(1:11, GO_size=BG_in, p.adjust_across=p.adjust_mine, p.adjust_MCls)    
   opfn <- paste(outdir, "GO_list/2_total_test_", MCl,"_", direction, ".xlsx", sep="")
   openxlsx::write.xlsx(x2, file=opfn, overwrite=T)

   ### significance GOs   
   x2 <- x%>%filter(Cluster==ii, BG_in>=5&BG_in<=500, (p.adjust<0.1)|(p.adjust_mine<0.1)|(p.adjust_MCls<0.1))
   x2 <- x2%>%arrange(pvalue)%>%
       dplyr::select(1:11, GO_size=BG_in, p.adjust_across=p.adjust_mine, p.adjust_MCls)
   ### 
   opfn <- paste(outdir, "GO_list/3_sigGO_", MCl, "_", direction, ".xlsx", sep="")
   openxlsx::write.xlsx(x2, file=opfn, overwrite=T) 
}    

### source data for figure 3 plots

fn <- "./5_enriched.outs/Correct_results/enriched_GO/combine_enriched_GO.xlsx"
x <- read.xlsx(fn)
opfn <- "./5_enriched.outs/Correct_results/TableS_Figure3_treeplot.xlsx"
write.xlsx(x, opfn)
    



######################################
### compare correct vs old results ###
######################################

## cg <- read_rds("./5_enriched.outs/Correct_results/1_enriched.direction.rds")
## cg2 <- cg%>%mutate(maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))%>%
##     filter(maxGSSize>=5&maxGSSize<=500)
## x <- cg2@compareClusterResult

## ###
## old <- read_rds("./5_enriched.outs/Old_results/1_enriched.direction.rds")
## old2 <- old%>%mutate(maxGSSize=as.numeric(gsub("/.*", "", BgRatio)))%>%
##     filter(maxGSSize>=5&maxGSSize<=500)
## x_old <- old2@compareClusterResult

## ###
## conditions <- as.character(sort(unique(x$Cluster)))

## summDF <- map_dfr(conditions, function(ii){
##     ##    
##     cat(ii, "\n")
##     ###
##     x2 <- x%>%dplyr::filter(Cluster==ii)%>%mutate(is_GO=ifelse(p.adjust<0.1, 1, 0))
##     GO <- x2%>%filter(is_GO==1)%>%pull(ID)
    
##     x_old2 <- x_old%>%dplyr::filter(Cluster==ii)%>%mutate(is_GO=ifelse(p.adjust<0.1, 1, 0))
##     GO_old <- x_old2%>%filter(is_GO==1)%>%pull(ID)
    
##     ##
##     summ <- data.frame(Cluster=ii,
##     Cluster2=ifelse(grepl("_1", ii), gsub("_1", "_Up", ii), gsub("_2", "_Down", ii)),
##     ntotal=nrow(x2), nsig=sum(x2$is_GO),
##     ntotal_old=nrow(x_old2), nsig_old=sum(x_old2$is_GO),
##     nolap=length(intersect(GO, GO_old)))
##     summ
## })    

## ##
## opfn <-"./5_enriched.outs/0_compare.xlsx"
## openxlsx::write.xlsx(summDF, opfn, overwrite=T)





