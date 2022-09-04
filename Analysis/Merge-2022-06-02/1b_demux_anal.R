#
library(tidyverse)
library(parallel)
##library(data.table)

### 
outFolder="./1_demux2_output/"

opfn <- "./1_demux2_output/1_demux_New.SNG.rds" ## 150,969
demux <- read_rds(opfn)

aa <- demux %>% dplyr::filter(NUM.READS>20,NUM.SNPS>20,NUM.READS<20000) %>%
    select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,Sample_ID=SNG.BEST.GUESS)


aa %>% group_by(EXP,Sample_ID) %>%
    summarize(n=n()) %>%
    filter(n>100) %>%
    spread(EXP,n) %>%
    write_tsv(paste0(outFolder,"cell.counts.matrix.tsv"))



cell.counts <- aa %>% group_by(EXP,Sample_ID) %>%
    summarize(n=n()) %>%
    filter(n>200)
cell.counts

#spread(cell.counts,EXP,n)
write_tsv(cell.counts,paste0(outFolder,"cell.counts.tsv"))


## END HERE

library(Seurat)
## library(SeuratDisk)
## library(SeuratData)
##
library(ggrastr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())


outdir <- "./1_demux2_v2_output/"

demux <- read_rds("./1_demux2_v2_output/1_demux_New.ALL.rds")
demux <- demux%>%dplyr::select("DROPLET.TYPE", "SNG.BEST.GUESS", "NEW_BARCODE", "EXP")

sc <- read_rds("./2_seurat.outs/1.0_seurat.all.rds")
meta <- sc@meta.data

x <- demux%>%dplyr::filter(NEW_BARCODE%in%meta$NEW_BARCODE)



####################
### droplet type ###
####################


dd <- x%>%
   group_by(EXP, DROPLET.TYPE)%>%
   summarise(ncell=n(),.groups="drop")
p <- ggplot(dd)+
   geom_bar(stat="identity", position=position_fill(reverse=F),
        aes(x=EXP, y=ncell, fill=DROPLET.TYPE))+
   theme_bw()+
   theme(legend.title=element_blank(),
         axis.title=element_blank(),
         axis.text.x=element_text(angle=90, hjust=1, size=8))

figfn <- paste(outdir, "Figure1.1_droplet.png", sep="")
png(figfn, width=420, height=420, res=120)
p
dev.off()



#######################################
### heatmap show individuals in EXP ###
#######################################

dd <- x%>%group_by(SNG.BEST.GUESS, EXP)%>%summarise(ncell=n(), .groups="drop")
dd <- dd%>%mutate(sampleID=SNG.BEST.GUESS)
###
dd2 <- dd%>%group_by(SNG.BEST.GUESS)%>%slice(which.max(ncell))%>%
    mutate(comb=paste(EXP, sampleID, sep="_"))
dd2 <- dd2%>%as.data.frame()%>%dplyr::select(sampleID, comb)

dd <- dd%>%left_join(dd2, by="sampleID")
 


### data for heatmap
dd3 <- dd%>%group_by(SNG.BEST.GUESS)%>%
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
dd4 <- dd%>%group_by(EXP)%>%
   mutate(Perc=ncell/sum(ncell))%>%as.data.frame()

### x labels
## sampleID <- dd2$sampleID
## names(sampleID) <- dd2$comb
p <- ggplot(dd4, aes(x=comb, y=EXP, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_blank(), #element_text(hjust=1, vjust=0.5, angle=90, size=6),
         axis.text.y=element_text(size=8),
         axis.title=element_blank())

###
figfn <- paste(outdir, "Figure1.3_heatmap.png", sep="")
png(figfn, width=800, height=320, res=100)
print(p)
dev.off()



####
####
dd0 <- dd%>%arrange(comb)
dd5 <- dd0[,1:3]%>%pivot_wider(names_from=EXP, values_from=ncell)
opfn <- paste(outdir, "1_EXP.sampleID.tsv", sep="")
write_tsv(dd5, opfn)






##SNG.BEST.GUESS
## cc <- read_tsv("../Covid19.Samples.txt") %>%
##     mutate(Sample_ID=paste(Pregnancy_ID,Origin,sep="-"))
## head(cc)

## aa <- aa %>% select(EXP,Sample_ID) %>% left_join(cc) %>%
##     select(EXP,Sample_ID,Pregnancy_ID,Origin,Condition)#,FetalSex)


## cell.counts <- aa %>% group_by(EXP,Sample_ID,Pregnancy_ID,Origin,Condition) %>%
##     summarize(n=n()) %>%
##     filter(n>100)

## write_tsv(cell.counts,paste0(outFolder,"cell.counts.tsv"))



#####cc <- read_tsv("../Covid19.Samples.v1.txt")
# which(cc$Pregnancy_ID %in% c("HPL20919","HPL20921","HPL20922","HPL20928","HPL20929"))

