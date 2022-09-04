######################################
### demux2 for downstream analysis ###
######################################
library(Matrix)
library(tidyverse)
library(parallel)
library(data.table)



rm(list=ls())

outdir <- "./1_demux_output/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)


############################
### 2, read demuxlet data ###
############################

# new column "NEW_BARCODE" was added to the data frame associated to each experiment (experiment name+barcode) "-1" was removed from end of barcode

basefolder <- "../../counts_brain_2022_merge/demuxlet/demuxlet/"
demuxfn <- list.files(basefolder, pattern="*.out.best")

# all the libraries from the membrane placenta 
expNames <- gsub(".out.best", "", demuxfn) 
expNames

## remove bad emulsions
demux <- mclapply(expNames,function(ii){
  cat("#Loading ", ii, "\n")
  fn <- paste0(basefolder, ii,".out.best")
  dd <- data.frame(fread(fn,header=T))
  dd <- dd%>%mutate(NEW_BARCODE=paste0(ii,"_", gsub("-1", "", BARCODE)),EXP=ii)
},mc.cores=6)


# merging all the experiments  (12 experiments)
demux <- do.call(rbind,demux)
###
dim(demux)

##



## saving the final demux file obtained from all experiments into "1_demux_New.ALL.rds" file
### output
opfn <- paste(outdir, "1_demux_New.ALL.rds", sep="")
write_rds(demux, opfn)

## DROPLET.TYPE: "SNG" 
opfn <- paste(outdir, "1_demux_New.SNG.rds", sep="") 
demux %>% filter(DROPLET.TYPE=="SNG") %>%
    write_rds(opfn)


###############################
### summary demulet results ###
###############################

demux <- read_rds("./1_demux_output/1_demux_New.SNG.rds")
aa <- demux %>% dplyr::filter(NUM.READS>20,NUM.SNPS>20,NUM.READS<20000) %>%
    select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,Sample_ID=SNG.BEST.GUESS)

##
aa2 <- aa %>% group_by(EXP,Sample_ID) %>%
    summarize(n=n(), .groups="drop")%>%ungroup()%>%
    pivot_wider(id_cols=EXP, names_from=Sample_ID, values_from=n, values_fill=0)

###
opfn <- paste(outdir, "cell.counts.matrix.tsv", sep="")
write_tsv(aa2, opfn)



aa3 <- aa %>% group_by(EXP,Sample_ID) %>%
    summarize(n=n(), .groups="drop")%>%ungroup()
opfn <- paste(outdir, "cell.counts.tsv", sep="")
write_tsv(aa3, opfn)




##
##

### droplet type
## dd <- x%>%
##    group_by(EXP, DROPLET.TYPE)%>%
##    summarise(ncell=n(),.groups="drop")
## p <- ggplot(dd)+
##    geom_bar(stat="identity", position=position_fill(reverse=F),
##         aes(x=EXP, y=ncell, fill=DROPLET.TYPE))+
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          axis.title=element_blank(),
##          axis.text.x=element_text(angle=90, hjust=1, size=8))

## figfn <- paste(outdir, "Figure1.1_droplet.png", sep="")
## png(figfn, width=420, height=420, res=120)
## p
## dev.off()



#######################################
### heatmap show individuals in EXP ###
#######################################

###
### heatmap 
demux <- read_rds("./1_demux_output/1_demux_New.SNG.rds")
demux2 <- demux %>% dplyr::filter(NUM.READS>20,NUM.SNPS>20,NUM.READS<20000) %>%
    select(NEW_BARCODE, NUM.READS, NUM.SNPS, EXP, sampleID=SNG.BEST.GUESS)

dd <- demux2%>%group_by(sampleID, EXP)%>%
   summarise(ncell=n(), .groups="drop")%>%
   filter(ncell>100)

###
dd2 <- dd%>%group_by(sampleID)%>%top_n(1,wt=ncell)%>%
    mutate(comb=paste(EXP, sampleID, sep="_"))
dd2 <- dd2%>%as.data.frame()%>%dplyr::select(sampleID, comb)

dd <- dd%>%left_join(dd2, by="sampleID")
 


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
## dd0 <- dd%>%arrange(comb)
## dd5 <- dd0[,1:3]%>%pivot_wider(names_from=EXP, values_from=ncell)
## opfn <- paste(outdir, "1_EXP.sampleID.tsv", sep="")
## write_tsv(dd5, opfn)








#
