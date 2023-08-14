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
##
library(ggrastr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(viridis)


#############################
### plots for publication ###
#############################

## plots_pub2, reviews of paper
###
outdir <- "./3_cluster.outs/plots_pub2/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


sc <- read_rds("./3_cluster.outs/1_seurat.cluster14.rds")

MCls_name <- c("0"="ODC", "5"="ODC", "10"="ODC", "13"="ODC",
    "1"="Astrocyte", "2"="Microglia", "11"="Microglia", "3"="OPC",
     "4"="Non_DA", "6"="Pericyte", "7"="DA", "8"="Endothelial", "9"="T-cell", "12"="Ependymal")


meta2 <- sc@meta.data%>%mutate(celltype=MCls_name[as.character(seurat_clusters)])

sc <- AddMetaData(sc, meta2)

opfn <- "./3_cluster.outs/1.2_seurat.annot.rds"
write_rds(sc, opfn)



###
### umap

sc <- read_rds("./3_cluster.outs/1.2_seurat.annot.rds")

### color version 1, "#b3de69", "#fdb462",  "#bc80bd", 
col_MCls <- c("ODC"="#fb8072", "Astrocyte"="#92cd2d", "Microglia"="#c38dc4", "OPC"="#35978f",
   "Non_DA"="#fc9016",  "Pericyte"="#bebada", "DA"="#80b1d3", "Endothelial"="#fccde5",
    "T-cell"="#8dd3c7",  "Ependymal"="#828282")

### color version 2
## col_MCls <- c("ODC"="#4daf4a", "Astrocyte"="#377eb8", "Microglia"="#e41a1c",  "OPC"="#984ea3",
##    "Non_DA"="#ff7f00",  "Pericyte"="#f781bf", "DA"="#ffff33", "Endothelial"="#a65628",
##     "T-cell"="#35978f", "Ependymal"="#828282")
 
p0 <- DimPlot(sc, reduction="umap", group.by="celltype", label=F, label.size=3.5, repel=T, raster=FALSE)+
   theme_bw()+
   scale_color_manual(values=col_MCls)+
   xlab("UMAP_1")+ylab("UMAP_2")+
   theme(legend.position="none",
         plot.title=element_blank())

###
figfn <- paste(outdir, "Figure1.2.1_umap.MCls.png", sep="")
png(figfn, width=500, height=450, res=120)
print(p0)
dev.off()



#################
### dot plots ###
#################

sc <- read_rds("./3_cluster.outs/1.2_seurat.annot.rds")

geneList <- c("MOBP", "MBP", "PLP1", "CNP",
   "VCAN", "PDGFRA", "OLIG1", "OLIG2",            
   "GFAP", "AQP4", "SLC1A2", "SLC4A4",
   "CSF1R", "CX3CR1", "LRRK1", "DOCK8", "C3", "P2RY12",
   "SLC18A2", "DDC", "TH", "SLC6A3",    
   "GAD1", "GAD2", "SLC17A7",
   "NR4A2", "PDGFRB", "MYO1B", "RGS5",    
   "PECAM1", "KDR", "CDH5", "FLT1", "CLDN5", "PTPRB", 
   "SKAP1", "THEMIS", "IL7R", "CD96",
   "HTR2C", "TTR")
## "CSPG4",  
gene_val <- c(rep(1,4), rep(2, 4), rep(3,4), rep(4,6), rep(5,4), rep(6, 3),
    rep(7,4), rep(8, 6), rep(9,4), rep(10, 2))
names(gene_val) <- geneList


## MCls <- c("0"="0:ODC", "5"="0:ODC", "10"="0:ODC", "13"="0:ODC",
##           "1"="1:Astrocyte", "2"="2:Microglia", "11"="2:Microglia",
##           "3"="3:OPC",
##           "4"="4:Non_DA",
##           "7"="5:DA",
##           "6"="6:Pericyte", "8"="7:Endothelial",
##           "9"="8:T-cell", "12"="9:Ependymal")

## sc2 <- RenameIdents(sc, MCls)


## MCls_val <- c("ODC"=1, "Astrocyte"=2, "Microglia"=3, "OPC"=4, "Non_DA"="5", "DA"=6,
##     "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)

## p <- DotPlot(sc2, features=geneList, cols=c("lightgrey", "red"),
##    scale=F, dot.scale=8)+
##     theme(legend.title=element_text(size=9),
##           axis.title=element_blank(),
##           axis.text.y=element_text(size=9),
##           axis.text.x=element_text(angle=90, hjust=1, size=8))
## ###
## figfn <- paste(outdir, "Figure2.1_dotplot.png", sep="")
## png(figfn, width=950, height=500, res=120)
## print(p)
## dev.off()

##

MCls_val <- c("ODC"=1, "OPC"=2, "Astrocyte"=3, "Microglia"=4, "DA"=5, "Non_DA"=6, 
    "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)

data <- FetchData(sc, vars=geneList)

meta <- sc@meta.data 
MCls <- sort(unique(meta$celltype))
 
###
DF_plot <- map_dfr(MCls, function(ii){
   ##
   cellSel <- meta%>%filter(celltype==ii)%>%pull(NEW_BARCODE)
   xi <- data[cellSel,]
   DF <- data.frame(MCls=ii, gene=colnames(data),
      yave=apply(xi, 2, mean), percent=apply(xi>0, 2, mean))
   ##DF <- DF%>%mutate(yave2=ifelse(percent<1e-6, yave, yave/percent))
   DF 
})
rownames(DF_plot) <- NULL

DF_plot <- DF_plot%>%group_by(gene)%>%
    mutate(yave2=yave/max(yave), percent2=percent/max(percent))%>%as.data.frame()


DF_plot2 <- DF_plot%>%mutate(MCls_value=as.numeric(MCls_val[MCls]), gene_value=as.numeric(gene_val[gene]))
DF_plot2 <- DF_plot2%>%
   mutate(MCl2=fct_reorder(MCls, MCls_value),
          gene2=fct_reorder(gene, gene_value),
          size2=ifelse(percent2<0.6, NA, percent2),
          yave2=ifelse(yave2<0.6, NA, yave2),
          gr2=ifelse(is.na(yave2), NA, "gr2"))
  
p <- ggplot(DF_plot2, aes(x=gene2, y=MCl2))+
     geom_point(aes(fill=yave2, color=gr2), size=4.5,  shape=21)+
     scale_fill_gradient(name="Average expression",
          low="white", high="#de2d26", na.value=NA,
          breaks=waiver(), n.breaks=5, limits=c(0,1),
          guide=guide_colourbar(barwidth=grid::unit(0.8,"lines"),
          barheight=grid::unit(5,"lines")))+
##    scale_size_binned(range=c(0.1, 7), limits=c(0.6,1), guide="none")+
    scale_color_manual(values=c("gr2"="black"), na.value=NA,  guide="none")+    
     ##scale_size_manual(values=c("a1"=0.1, "a2"=4), guide="none")+ 
    theme(axis.title=element_blank(),
          ##axis.text.y=element_text(size=9),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          panel.background=element_blank(),
          panel.border=element_rect(color="black", fill=NA, size=1.5),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          ##axis.line=element_line(color="black", size=2),
          legend.title=element_blank())
          #legend.text=element_text(size=8))

##
figfn <- paste(outdir, "Figure2.1_dotplot.png", sep="")
png(figfn, width=1000, height=500, res=120)
print(p)
dev.off()


### output data for plots
dot_df <- DF_plot2%>%
    dplyr::select(gene, gene_value, MCls, MCls_value, yave, yave2, percent=percent2)

opfn <- paste(outdir, "TableS_Figure1_dotplot.txt", sep="")
write.table(dot_df, file=opfn, row.names=F, col.names=T, quote=F)



#####################
### Violin plots ####
#####################

## col_MCls <- c("ODC"="#4daf4a", "Astrocyte"="#377eb8", "Microglia"="#e41a1c",  "OPC"="#984ea3",
##    "Non_DA"="#ff7f00",  "Pericyte"="#f781bf", "DA"="#ffff33", "Endothelial"="#a65628",
##     "T-cell"="#35978f", "Ependymal"="#828282")

## geneList <- c("MOBP", "MBP", 
##     "AQP4", "SLC1A2",
##    "CSF1R", "CX3CR1",  
##    "VCAN", "PDGFRA", 
##    "GAD1", "GAD2", 
##    "SLC18A2", "TH", 
##    "MYO1B", "RGS5",    
##    "FLT1", "CLDN5", 
##    "IL7R", "CD96",
##    "HTR2C", "TTR")

## cols <- c("#377eb8", "#ffff33", "#a65628", "#828282",
##           "#e41a1c", "#ff7f00", "#4daf4a", "#984ea3", "#f781bf", "#35978f")

## MCls_val <- c("ODC"=1, "Astrocyte"=2, "Microglia"=3, "OPC"=4, "Non_DA"="5", "DA"=6,
##     "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)

## data <- FetchData(sc, vars=geneList)

## meta <- sc@meta.data 
## MCls <- sort(unique(meta$celltype))
 
## ###
## DF_plot <- map_dfr(MCls, function(ii){
##    ##
##    cellSel <- meta%>%filter(celltype==ii)%>%pull(NEW_BARCODE)
##    xi <- data[cellSel,]
##    DF <- data.frame(MCls=ii, gene=colnames(data),
##       yave=apply(xi, 2, mean), percent=apply(xi>0, 2, mean))
##    ##DF <- DF%>%mutate(yave2=ifelse(percent<1e-6, yave, yave/percent))
##    DF 
## })
## rownames(DF_plot) <- NULL

 
## ##
## figfn <- paste(outdir, "Figure2.2_violin.png", sep="")
## png(figfn, width=800, height=550)
## print(vln)
## dev.off()



###
###

## meta <- sc@meta.data
## meta <- meta%>%mutate(Batch=ifelse(grepl("AKB", sampleID), "Bannon", "Mash"))

## x <- read.csv("BrainCV_cell.counts.final.csv")%>%filter(USE==1)%>%
##     dplyr::select(sampleID, Category)

## meta <- meta%>%inner_join(x, by="sampleID")

## MCls_val <- c("ODC"=1, "Astrocyte"=2, "Microglia"=3, "OPC"=4, "Non_DA"="5", "DA"=6,
##     "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)

## ###
## plot_DF <- meta%>%group_by(celltype, Category)%>%summarise(ncell=n(),.groups="drop")
## plot_DF <- plot_DF%>%group_by(celltype)%>%mutate(percent=ncell/sum(ncell))%>%ungroup()
## plot_DF2 <- plot_DF%>%mutate(MCls_value=as.numeric(MCls_val[celltype]),
##     MCl2=fct_reorder(celltype, MCls_value), comb=paste(celltype, Category, sep="_"))                       


## ###
## anno_df <- meta%>%group_by(celltype, Category)%>%
##    distinct(sampleID, .keep_all=T)%>%summarize(nindi=n(),.groups="drop")%>%as.data.frame()%>%
##    mutate(comb=paste(celltype, Category, sep="_"))%>%dplyr::select(nindi, comb) 
   
## plot_DF2 <- plot_DF2%>%left_join(anno_df, by="comb")

 
## p <- ggplot(plot_DF2, aes(x=MCl2, y=percent, fill=Category, group=Category))+
##    geom_bar(stat="identity", position="dodge")+
##    geom_text(aes(label=nindi), vjust=1.2, size=2.5)+ 
##    scale_fill_manual(values=c("Control"="#4575b4", "Opioid"="#d73027"))+
##    ylab("Percentage")+ 
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          axis.title.x=element_blank(),
##          axis.text.x=element_text(angle=90, hjust=0.5, size=10))

## figfn <- paste(outdir, "Figure3.1_barplot.png", sep="")
## png(figfn, width=880, height=420, res=120)
## print(p)
## dev.off()



##############################################
#### barplot show percentage of cell-types ###
##############################################

sc <- read_rds("./3_cluster.outs/1.2_seurat.annot.rds")


meta <- sc@meta.data


### cell-type proportion  
## summ <- meta%>%group_by(celltype)%>%summarise(ncell=n(), nt=nrow(meta), prop=ncell/nt, .groups="drop")%>%ungroup()
## summ <- summ%>%mutate(MCls_value=as.numeric(MCls_val[celltype]))%>%arrange(MCls_value)
## summ2 <- summ%>%dplyr::select(celltype, ncell, proportion=prop)
## ###
## opfn <- paste(outdir,  "Table1_celltype_prop.xlsx", sep="")
## openxlsx::write.xlsx(summ2, file=opfn, overwrite=T)


## x <- read.csv("BrainCV_cell.counts.final.csv")%>%filter(USE==1)%>%
##     dplyr::select(sampleID, Category) ## obtain treatment 
  
x <- openxlsx::read.xlsx("BrainCV_cell_counts_final_correct_2022-12-06.xlsx")%>%
   dplyr::select(sampleID, Category)


df2 <- meta%>%group_by(celltype, sampleID)%>%summarise(ncell=n(),.groups="drop")%>%ungroup()
df2 <- df2%>%group_by(sampleID)%>%mutate(total=sum(ncell), percet=ncell/total)%>%ungroup()
df2 <- df2%>%left_join(x, by="sampleID")

opfn <- paste(outdir, "Fig1_4_plotdata.rds", sep="")
write_rds(df2, opfn)

fn <- paste(outdir, "Fig1_4_plotdata.rds", sep="")
df2 <- read_rds(fn)


###
plotDF <- df2%>%group_by(celltype, Category)%>%
   summarise(ymean=mean(percet,na.rm=T), sd=sd(percet, na.rm=T), .groups="drop")


MCls_val <- c("ODC"=1, "OPC"=2, "Astrocyte"=3, "Microglia"=4, "DA"=5, "Non_DA"=6, 
    "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)
## col_MCls <- c("ODC"="#fb8072", "Astrocyte"="#b3de69", "Microglia"="#bebada", "OPC"="#8dd3c7",
##    "Non_DA"="#fdb462",  "Pericyte"="#ffffb3", "DA"="#80b1d3", "Endothelial"="#fccde5",
##     "T-cell"="#35978f", "Ependymal"="#828282")


plotDF <- plotDF%>%
    mutate(ylower=ymean-sd, ylower=ifelse(ylower<0, 0, ylower),
           yupper=ymean+sd,
           MCl_value=as.numeric(MCls_val[as.character(celltype)]),
           celltype2=fct_reorder(celltype, MCl_value, .desc=F))


p <- ggplot(plotDF, aes(x=celltype2, y=ymean, fill=Category, group=Category))+
   geom_bar(stat="identity", position=position_dodge(), alpha=0.8)+
   geom_errorbar(aes(ymin=ylower, ymax=yupper, color=Category), width=0.3, size=0.5, position=position_dodge(0.9))+ 
   scale_fill_manual(values=c("Control"="#0571b0", "Opioid"="#ca0020"))+
   scale_color_manual(values=c("Control"="#0571b0", "Opioid"="#ca0020"), guide="none")+ 
   ## scale_alpha_manual(values=c("Control"=1, "Opioid"=0.6),
   ##    guide=guide_legend(override.aes=list(fill="#fb8072")))+   
   ylab("Propotion by cell type")+
   theme_bw()+
   theme(legend.position=c(0.7,0.8),
         legend.title=element_blank(),
         legend.key.size=grid::unit(0.8, "lines"),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.text=element_text(size=10),
         axis.title.x=element_blank(),
         axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10),         
         axis.title.y=element_text(size=12),
         axis.text.y=element_text(size=10))
 
figfn <- paste(outdir, "Fig1_4_barplot.png", sep="")
png(figfn, width=540, height=420, res=120)
print(p)
dev.off()

###
###
x <- plotDF%>%pivot_wider(id_cols=celltype, names_from=Category, values_from=c(ymean, sd))%>%
   mutate(MCls_value=as.numeric(MCls_val[celltype]))%>%arrange(MCls_value)

### t-test
celltype <- sort(unique(df2$celltype))
t_test <- map_dfr(celltype, function(ii){
   ###    
   x1 <- df2%>%filter(celltype==ii, Category=="Control")%>%pull(percet)
   x2 <- df2%>%filter(celltype==ii, Category=="Opioid")%>%pull(percet)
   t0 <- t.test(x1, x2)
   pval <- t0$p.value 
   if ( pval>0.05){
       symbol <- "NS"
   }else if (pval>0.01 &pval<0.05){
       symbol <- "*"
   }else if (pval>0.001&pval<0.01){
       symbol <- "**"
   }else{
       symbol <- "***"
   }    
   tmp <- data.frame(celltype=ii, pvalue=pval, symbol=symbol)
   tmp
})
x2 <- x%>%left_join(t_test, by="celltype")

### output
opfn <- paste(outdir, "Table2_celltype_prop_OvsC.xlsx", sep="")
openxlsx::write.xlsx(x2%>%dplyr::select(-MCls_value), file=opfn, overwrite=T)


##################################
### percentage of boxplot plots ###
##################################

###
MCls_val <- c("ODC"=1, "OPC"=2, "Astrocyte"=3, "Microglia"=4, "DA"=5, "Non_DA"=6, 
    "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)

xlab_new <- c("ODC"="ODC", "OPC"="OPC", "Astrocyte"="Astrocyte", "Microglia"="Microglia",
   "DA"="DA neuron", "Non_DA"="Non DA neuron",
   "Pericyte"="Pericyte", "Endothelial"="Endothelia", "T-cell"="T-cell", "Ependymal"="Ependymal")       

###
fn <- paste(outdir, "Fig1_4_plotdata.rds", sep="")
df2 <- read_rds(fn)

df2 <- df2%>%dplyr::rename("percent"="percet")
### 
df2 <- df2%>%mutate(MCl_value=as.numeric(MCls_val[as.character(celltype)]),
                    celltype2=fct_reorder(celltype, MCl_value))
 
###  
p2 <- ggplot(df2, aes(x=celltype2, y=percent))+
   ## geom_violin(aes(fill=Category), alpha=0.8, width=1.2, position="dodge")+
   geom_boxplot(aes(fill=Category), lwd=0.25, position=position_dodge(0.9),
                outlier.fill=NA, outlier.shape=21, outlier.size=1)+
   ## stat_boxplot(geom="errorbar", width=0.25, position=position_dodge(0.75))+
   ## stat_summary(fun=median, color="grey", position=position_dodge(0.75),
   ##              geom="point", shape=23, size=3, show.legend=F)+
   scale_x_discrete(labels=xlab_new)+ 
   ylab("Propotion by cell type")+
   scale_fill_manual(values=c("Control"="#0571b0", "Opioid"="#ca0020"))+
   ## scale_color_manual(values=c("Control"="#0571b0", "Opioid"="#ca0020"), guide="none")+   
   theme_bw()+
   theme(legend.position=c(0.7, 0.8),
         legend.title=element_blank(),
         legend.text=element_text(size=10),
         ##legend.key.size=grid::unit(0.6, "lines"),
         legend.background=element_blank(),
         legend.box.background=element_blank(),         
         axis.title.x=element_blank(),
         axis.title.y=element_text(size=12),
         axis.text.x=element_text(angle=45, hjust=1, size=10),
         axis.text.y=element_text(size=10))   


figfn <- paste(outdir, "Fig1_4.2_boxplot.png", sep="")
png(figfn, width=540, height=420, res=120)
print(p2)
dev.off()

### output data
fn <- "./3_cluster.outs/plots_pub2/Fig1_4_plotdata.rds"
df2 <- read_rds(fn)
df2 <- df2%>%dplyr::rename(percent=percet)

opfn <- paste(outdir, "TableS_Figure1_boxplot.txt", sep="")
write.table(df2, file=opfn, row.names=F, col.names=T, quote=F)
 



## ###
## ###
## p3 <- ggplot(df2, aes(x=celltype2, y=log2(percent)))+
##    geom_violin(aes(fill=Category), alpha=0.8, position="dodge")+
##    ##geom_boxplot(aes(fill=Category), width=0.8, outlier.shape=NA)+
##    ylab("Percents of cell-types in each sample")+
##    scale_fill_manual(values=c("Control"="#0571b0", "Opioid"="#ca0020"))+
##    ##scale_color_manual(values=c("Control"="#0571b0", "Opioid"="#ca0020"), guide="none")+ 
##    theme_bw()+
##    theme(legend.position=c(0.8,0.8),
##          legend.title=element_blank(),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          axis.title.x=element_blank(),
##          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10),
##          axis.text.y=element_text(size=10))   

## figfn <- paste(outdir, "Figure3.4_violin.png", sep="")
## png(figfn, width=580, height=420, res=120)
## print(p3)
## dev.off()



##################################################
### umap plots facet by cohorts and treatments ###
##################################################

sc <- read_rds("./3_cluster.outs/1.2_seurat.annot.rds")

###
###
## MCls_val <- c("ODC"=1, "Astrocyte"=2, "Microglia"=3, "OPC"=4, "Non_DA"="5", "DA"=6,
##     "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)
 
MCls_val <- c("ODC"=1, "OPC"=2, "Astrocyte"=3, "Microglia"=4, "DA"=5, "Non_DA"=6, 
    "Pericyte"=7, "Endothelial"=8, "T-cell"=9, "Ependymal"=10)

col_MCls <- c("ODC"="#fb8072", "OPC"="#35978f", "Astrocyte"="#92cd2d", "Microglia"="#c38dc4", 
   "DA"="#80b1d3", "Non_DA"="#fc9016",  "Pericyte"="#bebada",  "Endothelial"="#fccde5",
    "T-cell"="#8dd3c7",  "Ependymal"="#828282")

## col_MCls <- c("ODC"="#4daf4a", "Astrocyte"="#377eb8", "Microglia"="#e41a1c",  "OPC"="#984ea3",
##    "Non_DA"="#ff7f00",  "Pericyte"="#f781bf", "DA"="#ffff33", "Endothelial"="#a65628",
##     "T-cell"="#35978f", "Ependymal"="#828282")

## sc2 <- read_rds("./3_cluster.outs/1_seurat.cluster.rds")
umap <- as.data.frame(sc[["umap"]]@cell.embeddings)%>%rownames_to_column(var="NEW_BARCODE")

meta <- sc@meta.data
meta <- meta%>%mutate(Batch=ifelse(grepl("AKB", sampleID), "Detroit", "Miami"))

## x <- read.csv("BrainCV_cell.counts.final.csv")%>%filter(USE==1)%>%
##     dplyr::select(sampleID, Category)

x <- openxlsx::read.xlsx("BrainCV_cell_counts_final_correct_2022-12-06.xlsx")%>%
   dplyr::select(sampleID, Category) 

meta <- meta%>%inner_join(x, by="sampleID")

 
df2 <- meta%>%dplyr::select(Batch, Category, NEW_BARCODE, celltype)%>%inner_join(umap, by="NEW_BARCODE")

## ###
p <- ggplot(df2, aes(x=UMAP_1, y=UMAP_2, colour=factor(celltype)))+
    geom_point(size=0.1)+
   ## rasterise(geom_point(size=0.1),dpi=300)+
   scale_colour_manual(values=col_MCls,
       guide=guide_legend(override.aes=list(size=2, ncol=1)))+
   facet_grid(Batch~Category, scales="fixed")+ 
   ## guides(col=guide_legend(override.aes=list(size=2), ncol=1))+
   theme_bw()+
   theme(legend.position="none",
         legend.title=element_blank(),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(1,"lines"),
         legend.text=element_text(size=10),
         strip.text=element_text(size=14),
         axis.text=element_text(size=10),
         axis.title=element_text(size=12))
###
figfn <- paste(outdir, "FigureS1.1_cluster.png", sep="")
png(figfn, width=720, height=720, res=120)
p
dev.off()

###
### output data for plots

meta2 <- meta%>%dplyr::select(NEW_BARCODE, Cohort=Batch, Condition=Category, Celltype=celltype)
meta2 <- meta2%>%inner_join(umap, by="NEW_BARCODE")
opfn2 <- paste(outdir, "TableS_Figure1_umap.txt", sep="")
write.table(meta2, file=opfn2, row.names=F, col.names=T, quote=F)







### 0, previous version, raster=T
### 1, current version, raster=F

###
### 
