#
library(tidyverse)
library(parallel)
##library(data.table)

### 
outFolder="./1_demux_output/"

opfn <- "./1_demux_output/1_demux_New.SNG.rds" ## 150,969
demux <- read_rds(opfn)

aa <- demux %>% dplyr::filter(NUM.READS>40,NUM.SNPS>40) %>%
   select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,Sample_ID=SNG.BEST.GUESS) %>% 
   mutate(Sample_ID2=gsub(".*_","",Sample_ID))


aa %>% group_by(EXP,Sample_ID2) %>%
    summarize(n=n()) %>%
    filter(n>100) %>%
    spread(EXP,n) %>%
    write_tsv(paste0(outFolder,"cell.counts.matrix.v2.tsv"))


cv <- read_tsv("BrainCV.txt")



cell.counts <- aa %>% group_by(EXP,Sample_ID2) %>%
    summarize(n=n()) %>%
    filter(n>200)
cell.counts

#spread(cell.counts,EXP,n)
write_tsv(cell.counts,paste0(outFolder,"cell.counts.v2.tsv"))


cc2 <- cv %>% full_join(cell.counts)
dim(cc2)

write_tsv(cc2,paste0(outFolder,"cell.counts.cc2.tsv"))  

#######################################
### heatmap show individuals in EXP ###
#######################################

###
### heatmap 

dd <- aa %>%group_by(Sample_ID2, EXP)%>%
   summarise(ncell=n(), .groups="drop")%>%
   filter(ncell>100)

###
dd2 <- dd%>%group_by(Sample_ID2)%>%top_n(1,wt=ncell)%>%
   mutate(comb=paste(EXP, Sample_ID2, sep="_"))
dd2 <- dd2%>%as.data.frame()%>%dplyr::select(Sample_ID2, comb)

dd <- dd%>%left_join(dd2, by="sampleID")

dd -> cell.counts %>% mutate(comb=paste0(EXP,"_",Sample_ID2))

### data for heatmap
dd3 <- dd%>%group_by(sampleID)%>%
   mutate(Perc=ncell/sum(ncell))%>%as.data.frame()

### x labels
## sampleID <- dd2$sampleID
## names(sampleID) <- dd2$comb
p <- ggplot(dd3, aes(x=comb, y=EXP, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
                       low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_blank(), #element_text(hjust=1, vjust=0.5, angle=90, size=6),
         axis.text.y=element_text(size=8),
         axis.title=element_blank())

###
figfn <- paste(outdir, "Figure1.2_heatmap.png", sep="")
png(figfn, width=600, height=280, res=100)
print(p)
dev.off()



###
dd4 <- cc2 %>% 
   mutate(Sample_ID3=paste0(ExpectedPool,"-",Sample_ID2)) %>%
   group_by(EXP) %>%
   mutate(Perc=n/sum(n)) %>% 
   arrange(EXP)
    
##   as.data.frame()

### x labels
## sampleID <- dd2$sampleID
## names(sampleID) <- dd2$comb
p <- ggplot(dd4, aes(x=Sample_ID3, y=EXP, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
                       low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(#axis.text.x=element_blank(), 
      axis.text.x=element_text(hjust=1, vjust=1, angle=45, size=8),
         axis.text.y=element_text(size=8),
         axis.title=element_blank())

###
figfn <- paste(outFolder, "Figure1_heatmap_v2.pdf", sep="")
pdf(figfn, width=14, height=6)
print(p)
dev.off()


