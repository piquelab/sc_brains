##
library(Matrix)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library(qvalue)
library(openxlsx)
##
## library(DESeq2)
## library(annotables)
###library(biobroom)
##library(Seurat)
## library(cowplot)
## library(ggrepel)
##

rm(list=ls())

outdir <- "./1_TWAS.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

###
### Twas Traits
### Last modified Sept-25-2022, by Julong wei


#############################
### 1. Multiple tissues   ###
###    Q-Q plots          ###
 ############################

fn <- "smultixcan-mashr-pvalues.tsv.gz"
gene_pval <- fread(fn, data.table=F)
gene_pval <- gene_pval%>%distinct(gene_name, .keep_all=T)%>%column_to_rownames(var="gene_name")
colnames(gene_pval) <- gsub("-", "_", colnames(gene_pval))

traits <- colnames(gene_pval)

target_traits <- c("PGC_ADHD_EUR_2017",
   "20546_3_Substances_taken_for_depression_Medication_prescribed_to_you_for_at_least_two_weeks", 
   "1558_Alcohol_intake_frequency",
   "1498_Coffee_intake")

## data for QQ plot
plotDF <- map_dfr(target_traits, function(ii){
    ###
    cat(ii, "\n")
    plotDF2 <- data.frame(pval=gene_pval[,ii], traits=ii, trait2=ii)%>%arrange(pval)
    ngene <- nrow(plotDF2)
    plotDF2 <- plotDF2%>%mutate(observed=-log10(pval), expected=-log10(ppoints(ngene)))
    ###
    if (grepl("20546",ii)){
      plotDF2$trait2 <- "20546_Substances_for_depression"
    }    
    plotDF2
})

### qq plots
p  <- ggplot(plotDF, aes(x=expected, y=observed))+
   ggrastr::rasterise(geom_point(size=0.3, colour="grey30"),dpi=300)+
   facet_wrap(~trait2, ncol=2, scales="free")+
   geom_abline(colour="red")+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(strip.text=element_text(size=9))

 
figfn <- paste(outdir, "Figure1.1_qq_multiple.png", sep="")
png(figfn, width=580, height=580, res=120)
print(p)
dev.off()


#########################
### 2. single tissue  ###
#########################

fn <- "./spredixcan/spredixcan-Brain_Substantia_nigra.tsv.gz"
twas <- fread(fn, data.table=F)
twas <- twas%>%distinct(gene_name, .keep_all=T)%>%column_to_rownames(var="gene_name")

###
traits <- colnames(twas)
gene <- rownames(twas)

## x <- sort(traits[-1])
## dd <- data.frame(traits=x)
## opfn <- "traits_list.txt"
## write.table(dd, file=opfn, row.names=F, col.names=F, quote=F)
 

## data for QQ plot
plotDF <- map_dfr(target_traits, function(ii){
    ###
    cat(ii, "\n")
    zval <- twas[, ii]
    pval <- 2*pnorm(abs(zval), lower.tail=F)
    plotDF2 <- data.frame(pval=pval, traits=ii, trait2=ii)%>%arrange(pval)
    ngene <- nrow(plotDF2)
    plotDF2 <- plotDF2%>%mutate(observed=-log10(pval), expected=-log10(ppoints(ngene)))
    ###
    if (grepl("20546",ii)){
      plotDF2$trait2 <- "20546_Substances_for_depression"
    }    
    plotDF2
})

### qq plots
p2  <- ggplot(plotDF, aes(x=expected, y=observed))+
   ggrastr::rasterise(geom_point(size=0.3, colour="grey30"),dpi=300)+
   facet_wrap(~trait2, ncol=2, scales="free")+
   geom_abline(colour="red")+
   xlab(bquote("Expected"~ -log[10]~"("~italic(plain(P))~")"))+
   ylab(bquote("Observed"~-log[10]~"("~italic(plain(P))~")"))+
   theme_bw()+
   theme(strip.text=element_text(size=9))

 
figfn <- paste(outdir, "Figure1.2_qq_single.png", sep="")
png(figfn, width=580, height=580, res=120)
print(p2)
dev.off()



## ### p-value
## fn <- "../TWAS/smultixcan-mashr-pvalues.tsv.gz"
## gene_pval <- fread(fn, data.table=F)
## gene_pval <- gene_pval%>%distinct(gene_name, .keep_all=T)%>%column_to_rownames(var="gene_name")
## colnames(gene_pval) <- gsub("-", "_", colnames(gene_pval))



## ### rcp
## fn <- "../TWAS/fastenloc-torus-rcp.tsv.gz"
## gene_rcp <- fread(fn, data.table=F)
## gene_rcp <- gene_rcp%>%distinct(gene_id, .keep_all=T)%>%column_to_rownames(var="gene_id")
## colnames(gene_rcp) <- gsub("-", "_", colnames(gene_rcp))

## gene <- intersect(intersect(rownames(twas), rownames(gene_pval)), rownames(gene_rcp))



################################################################################# 
## 3. summary of TWAS for traits of interest, Single Tissue, substantia nigra ###
###          z-score from substantia_nigra tissue                             ###
#################################################################################

fn <- "./spredixcan/spredixcan-Brain_Substantia_nigra.tsv.gz"
twas <- fread(fn, data.table=F)
twas <- twas%>%distinct(gene_name, .keep_all=T)%>%column_to_rownames(var="gene_name")

###
traits <- colnames(twas)
gene <- rownames(twas)
geneAnno <- bitr(gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)


gene_comb <- NULL
###
### ADHD
target_trait <- traits[grepl("ADHD", traits)]

zval <- twas[, target_trait]
pval <- 2*pnorm(abs(zval), lower.tail=F)
FDR <- qvalue(pval)$qvalues

gene_df <- data.frame(gene=gene, zscore=zval, pval=pval, FDR=FDR, traits=target_trait, category="ADHD")

gene_df2 <- gene_df%>%filter(FDR<0.1)%>%left_join(geneAnno, by=c("gene"="ENSEMBL"))

gene_comb <- rbind(gene_comb, gene_df2)

## summ <- data.frame(traits=target_trait, ngene=nrow(gene_df2), category="ADHD")
## DF_comb <- rbind(DF_comb, summ)


## geneID <- bitr(gene_df2$gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)


###
### parkinson disease
target_traits <- traits[grepl("parkinson|Parkinson", traits)]
gene_df2 <- lapply(target_traits, function(ii){
    ###     
    zval <- twas[, ii]
    pval <- 2*pnorm(abs(zval), lower.tail=F)
    FDR <- qvalue(pval)$qvalues
    
    ###
    gene_df <- data.frame(gene=gene, zscore=zval, pval=pval, FDR=FDR, traits=ii, category="parkinson")
    gene_df2 <- gene_df%>%filter(FDR<0.1)%>%left_join(geneAnno, by=c("gene"="ENSEMBL"))
    ##
    cat(ii, nrow(gene_df2), "\n")
    gene_df2
})
gene_df2 <- do.call(rbind, gene_df2)
##
gene_comb <- rbind(gene_comb, gene_df2)

## summ <- gene_df2%>%group_by(traits)%>%summarise(ngene=n(), .groups="drop")%>%
##     mutate(category="parkinson")%>%as.data.frame()
## DF_comb <- rbind(DF_comb, summ)


## gene_df2 <- gene_df2%>%group_by(gene)%>%slice_max(order_by=abs(zscore),n=1)

## geneID <- bitr(gene_df2$gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)

## gene_df2 <- gene_df2%>%left_join(geneID, by=c("gene"="ENSEMBL"))%>%
##    dplyr::rename(zscore_TWAS=zscore)



##
### Addiction
target_traits <- traits[grepl("Addiction|addiction|addicted|Addicted|abuse|Abuse|substance|Substance", traits)]
gene_df2 <- lapply(target_traits, function(ii){
    ###     
    zval <- twas[, ii]
    pval <- 2*pnorm(abs(zval), lower.tail=F)
    FDR <- qvalue(pval)$qvalues
    
    ###
    gene_df <- data.frame(gene=gene, zscore=zval, pval=pval, FDR=FDR, traits=ii, category="addiction")
    gene_df2 <- gene_df%>%filter(FDR<0.1)%>%left_join(geneAnno, by=c("gene"="ENSEMBL"))
    ##
    cat(ii, nrow(gene_df2), "\n")
    gene_df2
})
gene_df2 <- do.call(rbind, gene_df2)

gene_comb <- rbind(gene_comb, gene_df2)

## summ <- gene_df2%>%group_by(traits)%>%summarise(ngene=n(), .groups="drop")%>%
##     mutate(category="addiction")%>%as.data.frame()

## DF_comb <- rbind(DF_comb, summ)

## gene_df2 <- gene_df2%>%group_by(gene)%>%slice_max(order_by=abs(zscore),n=1)

## geneID <- bitr(gene_df2$gene, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)

## gene_df2 <- gene_df2%>%left_join(geneID, by=c("gene"="ENSEMBL"))%>%
##    dplyr::rename(zscore_TWAS=zscore)



###
### alcohol

target_traits <- traits[grepl("alcohol|Alcohol|drink|Drink|Glass|wine|drunk|cup", traits)]
## target_traits <- target_traits[c(1:3, 6,11,15, 23:27)]
gene_df2 <- lapply(target_traits, function(ii){
    ###     
    zval <- twas[, ii]
    pval <- 2*pnorm(abs(zval), lower.tail=F)
    FDR <- qvalue(pval)$qvalues
    
    ###
    gene_df <- data.frame(gene=gene, zscore=zval, pval=pval, FDR=FDR, traits=ii, category="alcohol")
    gene_df2 <- gene_df%>%filter(FDR<0.1)%>%left_join(geneAnno, by=c("gene"="ENSEMBL"))
    ##
    cat(ii, nrow(gene_df2), "\n")
    gene_df2
})
gene_df2 <- do.call(rbind, gene_df2)

gene_comb <- rbind(gene_comb, gene_df2)


## summ <- gene_df2%>%group_by(traits)%>%summarise(ngene=n(), .groups="drop")%>%
##     mutate(category="alcohol")%>%as.data.frame()

## DF_comb <- rbind(DF_comb, summ)


###
### opioid
target_traits <- traits[grepl("opioid|Opioid|Oplate|morphine|Morphine", traits)]
gene_df2 <- lapply(target_traits, function(ii){
    ###     
    zval <- twas[, ii]
    pval <- 2*pnorm(abs(zval), lower.tail=F)
    FDR <- qvalue(pval)$qvalues
    
    ###
    gene_df <- data.frame(gene=gene, zscore=zval, pval=pval, FDR=FDR, traits=ii, category="opioid")
    gene_df2 <- gene_df%>%filter(FDR<0.1) %>%left_join(geneAnno, by=c("gene"="ENSEMBL"))
    ##
    cat(ii, nrow(gene_df2), "\n")
    gene_df2
})
gene_df2 <- do.call(rbind, gene_df2)

gene_comb <- rbind(gene_comb, gene_df2)


## summ <- gene_df2%>%group_by(traits)%>%summarise(ngene=n(), .groups="drop")%>%
##     mutate(category="opioid")%>%as.data.frame()

## DF_comb <- rbind(DF_comb, summ)


###
### caffeine
target_traits <- traits[grepl("Caffeine|caffeine|coffee|Coffee", traits)]
gene_df2 <- lapply(target_traits, function(ii){
    ###     
    zval <- twas[, ii]
    pval <- 2*pnorm(abs(zval), lower.tail=F)
    FDR <- qvalue(pval)$qvalues
    
    ###
    gene_df <- data.frame(gene=gene, zscore=zval, pval=pval, FDR=FDR, traits=ii, category="caffeine")
    gene_df2 <- gene_df%>%filter(FDR<0.1)%>%left_join(geneAnno, by=c("gene"="ENSEMBL"))
    ##
    cat(ii, nrow(gene_df2), "\n")
    gene_df2
})
gene_df2 <- do.call(rbind, gene_df2)

gene_comb <- rbind(gene_comb, gene_df2)

## summ <- gene_df2%>%group_by(traits)%>%summarise(ngene=n(), .groups="drop")%>%
##     mutate(category="caffeine")%>%as.data.frame()

## DF_comb <- rbind(DF_comb, summ)


###
### painkiller
target_traits <- traits[grepl("painkiller|Painkiller|Pain|pain|Benzodiazepine|benzodiazepine", traits)]

## target_traits <- traits[grepl("amphetamine|Amphetamine", traits)]

gene_df2 <- lapply(target_traits, function(ii){
    ###     
    zval <- twas[, ii]
    pval <- 2*pnorm(abs(zval), lower.tail=F)
    FDR <- qvalue(pval)$qvalues
    
    ###
    gene_df <- data.frame(gene=gene, zscore=zval, pval=pval, FDR=FDR, traits=ii, category="painkiller")
    gene_df2 <- gene_df%>%filter(FDR<0.1)%>%left_join(geneAnno, by=c("gene"="ENSEMBL"))
    ##
    cat(ii, nrow(gene_df2), "\n")
    gene_df2
})
gene_df2 <- do.call(rbind, gene_df2)

gene_comb <- rbind(gene_comb, gene_df2) 

##
opfn <- paste(outdir, "1_TWAS_traits.xlsx", sep="")
openxlsx::write.xlsx(gene_comb, opfn, overwrite=T)




###
### New traits, Tobacco traits

gene_comb <- NULL

target_traits <- traits[grepl("smoke|smoking|Cigarette|cigarette|Nicotine|nicotine|Tobacco|tobacco|Cigar|cigar|Chew|chew|Vaping|vaping|Vape|vape|Packet|packet", traits)]

gene_df2 <- lapply(target_traits, function(ii){
    ###     
    zval <- twas[, ii]
    pval <- 2*pnorm(abs(zval), lower.tail=F)
    FDR <- qvalue(pval)$qvalues
    
    ###
    gene_df <- data.frame(gene=gene, zscore=zval, pval=pval, FDR=FDR, traits=ii, category="smoke")
    gene_df2 <- gene_df%>%filter(FDR<0.1)%>%left_join(geneAnno, by=c("gene"="ENSEMBL"))
    ##
    cat(ii, nrow(gene_df2), "\n")
    gene_df2
})
gene_df2 <- do.call(rbind, gene_df2)

gene_comb <- rbind(gene_comb, gene_df2)


###
### New traits, Marijuana

target_traits <- traits[grepl("Marijuana|marijuana|Cannabis|cannabis|THC|Cannabinoid|cannabinoid|cannabinol|Cannabinol|mariguana|Mariguana", traits)]

gene_df2 <- lapply(target_traits, function(ii){
    ###     
    zval <- twas[, ii]
    pval <- 2*pnorm(abs(zval), lower.tail=F)
    FDR <- qvalue(pval)$qvalues
    
    ###
    gene_df <- data.frame(gene=gene, zscore=zval, pval=pval, FDR=FDR, traits=ii, category="marijuana")
    gene_df2 <- gene_df%>%filter(FDR<0.1)%>%left_join(geneAnno, by=c("gene"="ENSEMBL"))
    ##
    cat(ii, nrow(gene_df2), "\n")
    gene_df2
})
gene_df2 <- do.call(rbind, gene_df2)

gene_comb <- rbind(gene_comb, gene_df2)


###
comb2 <- openxlsx::read.xlsx("./1_TWAS.outs/1_TWAS_traits.xlsx")
gene_comb <- rbind(comb2, gene_comb)

###
opfn <- "./1_TWAS.outs/1_TWAS_traits.xlsx"
openxlsx::write.xlsx(gene_comb, opfn, overwrite=T)



########################
### Add other traits ###
########################



x <- openxlsx::read.xlsx("./1_TWAS.outs/traits.xlsx")


trait2 <- traits[grepl("schizophrenia|Schizophrenia", traits)]
DF2 <- data.frame(traits=unique(trait2),category="scz")
    
trait2 <- traits[grepl("bipolar|Bipolar", traits)]
DF3 <- data.frame(traits=unique(trait2),category="BIP")

trait2 <- traits[grepl("depression|Depression", traits)]
DF4 <- data.frame(traits=unique(trait2),category="Depression")

trait2 <- traits[grepl("height|Height|BMI", traits)]
DF5 <- data.frame(traits=unique(trait2),category="height")

trait2 <- traits[grepl("CAD|coronary|Coronary", traits)]
DF6 <- data.frame(traits=unique(trait2), category="CAD")

trait2 <- traits[grepl("BMI", traits)]
DF7 <- data.frame(traits=unique(trait2), category="BMI")

DF <- rbind(DF2, DF3, DF4, DF5, DF6, DF7)

###
### 
nsig_gene <- sapply(1:nrow(DF), function(i){
    ###
    ii <- DF$traits[i]
    category2 <- DF$category[i]
    
    zval <- twas[, ii]
    pval <- 2*pnorm(abs(zval), lower.tail=F)
    FDR <- qvalue(pval)$qvalues 
    sum(FDR<0.1,na.rm=T)
})

DF$nsig <- nsig_gene

DF2 <- DF%>%filter(nsig_gene>5)

###
###
opfn <- "./1_TWAS.outs/trait3_summ.xlsx"
openxlsx::write.xlsx(DF2, file=opfn, overwrite=T)


###
###
fn <- "./1_TWAS.outs/trait3_summ.xlsx"
DF <- openxlsx::read.xlsx(fn)
DF2 <- DF[9:19,]

gene_df2 <- map_dfr(1:nrow(DF2), function(i){
    ###

    ii <- DF2$traits[i]
    category2 <- DF2$category[i]
    
    zval <- twas[, ii]
    pval <- 2*pnorm(abs(zval), lower.tail=F)
    FDR <- qvalue(pval)$qvalues
    
    ###
    gene_df <- data.frame(gene=gene, zscore=zval, pval=pval, FDR=FDR, traits=ii, category=category2)
    gene_df2 <- gene_df%>%filter(FDR<0.1)%>%left_join(geneAnno, by=c("gene"="ENSEMBL"))
    ##
    cat(ii, nrow(gene_df2), "\n")
    gene_df2
})

## gene_comb <- openxlsx::read.xlsx("./1_TWAS.outs/1_TWAS_traits.xlsx")
## gene_comb <- rbind(gene_comb, gene_df2)

opfn <- "./1_TWAS.outs/1.2_TWAS_traits.xlsx"
openxlsx::write.xlsx(gene_df2, opfn, overwrite=T) 

##
###gene_comb <- rbind(gene_comb, gene_df2)




## summ <- gene_comb%>%group_by(traits, category)%>%summarise(ngene=n(), .groups="drop")

## opfn <- paste(outdir, "trait2_summ.xlsx", sep="")
## openxlsx::write.xlsx(summ, file=opfn, overwrite=T)



## summ <- gene_df2%>%group_by(traits)%>%summarise(ngene=n(), .groups="drop")%>%
##     mutate(category="pain")%>%as.data.frame()

## DF_comb <- rbind(DF_comb, summ)

## gene_comb2 <- gene_comb%>%left_join(resMat, by="SYMBOL")
## ###%>%left_join(trait_DF, by="traits")%>%filter(Relevant==1)
 
## summ2 <- gene_comb2%>%group_by(traits, category)%>%
##    summarise(ntwas=n(), olap_CL0=sum(Cluster0,na.rm=T), olap_CL1=sum(Cluster1,na.rm=T),
##       olap_CL2=sum(Cluster2,na.rm=T), olap_CL3=sum(Cluster3,na.rm=T), olap_CL4=sum(Cluster4,na.rm=T), olap_CL6=sum(Cluster6,na.rm=T), .groups="drop")

## summ2 <- summ2%>%arrange(category)

## ###
## ###
## opfn <- paste(outdir, "1.2_summary_TWAS_traits.xlsx", sep="")
## openxlsx::write.xlsx(summ2, opfn, overwrite=T)













