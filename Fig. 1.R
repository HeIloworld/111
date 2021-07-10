library(Seurat)
library(ggplot2)
library(monocle)
setwd("~/Desktop/R code/")

# =============================================================================================
#All sequenced cell, before removing blood cells
#==============================================================================================
expresswt2336 <- read.csv("wtexpress2336.csv",stringsAsFactors=F)
dim(expresswt2336)
rownames(expresswt2336) <- expresswt2336[,1]
rownames(expresswt2336)
grep("ERCC",rownames(expresswt2336),value = T)
expresswt2336 <- expresswt2336[setdiff(rownames(expresswt2336),grep("ERCC",rownames(expresswt2336),value = T)),]
expresswt2336 <- expresswt2336[,-1]

metadatawt2336 <- read.csv("wtmeta2336.csv")
rownames(metadatawt2336) <- metadatawt2336[,1]
colnames(metadatawt2336)
table(rownames(metadatawt2336)==colnames(expresswt2336))

wt2336 <- CreateSeuratObject(counts = expresswt2336, project = "wt2336", min.cells = 3, min.features = 1000) 
VlnPlot(wt2336, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(wt2336, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
wt2336 <- subset(wt2336, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & nCount_RNA < 200000)
wt2336 <- AddMetaData(wt2336,metadata = metadatawt2336 )
dim(wt2336@meta.data)

wt2336 <- NormalizeData(wt2336, normalization.method = "LogNormalize", scale.factor = 10000)
wt2336 <- FindVariableFeatures(wt2336, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(wt2336), 10)
plot1 <- VariableFeaturePlot(wt2336)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(wt2336)
wt2336 <- ScaleData(wt2336, features = all.genes) 
wt2336 <- RunPCA(wt2336, features = VariableFeatures(object = wt2336))
print(wt2336[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(wt2336, dims = 1:2, reduction = "pca")
DimPlot(wt2336, reduction = "pca")
DimHeatmap(wt2336, dims = 1, cells = 500, balanced = TRUE)

wt2336 <- FindNeighbors(wt2336, dims = 1:15)
wt2336 <- FindClusters(wt2336, resolution = 0.2)
wt2336 <- RunTSNE(wt2336, dims = 1:15,seed.use = 1,perplexity=50)
DimPlot(wt2336, reduction = "tsne",pt.size = 2,shape.by ="condition",label = T)

#gata1, blood cell
FeaturePlot(wt2336, features =c("gata1a","hbbe1.3","hbbe2"),cols = c("grey", "red"),pt.size=1,min.cutoff =0)

bloodcell <- WhichCells(wt2336,ident = 1)

#marker
FeaturePlot(wt2336, features =c("myh6","vmhc","fli1a","tcf21","twist1a","vim"),cols = c("grey", "red"),pt.size=2,min.cutoff =0)

# =============================================================================================
#Fig1b, remove blood cells, retain 1581 cells
#==============================================================================================
wt2336blood <- subset(wt2336, idents = 1,invert=T)
dim(wt2336blood@meta.data)
wt2336blood <- NormalizeData(wt2336blood, normalization.method = "LogNormalize", scale.factor = 10000)
wt2336blood <- FindVariableFeatures(wt2336blood, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(wt2336blood), 10)
plot1 <- VariableFeaturePlot(wt2336blood)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(wt2336blood)
wt2336blood <- ScaleData(wt2336blood, features = all.genes) #??Ҫvars.to.regress = "run" ??Ⱥ???ÿ?
wt2336blood <- RunPCA(wt2336blood, features = VariableFeatures(object = wt2336blood))
print(wt2336blood[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(wt2336blood, dims = 1:2, reduction = "pca")
DimPlot(wt2336blood, reduction = "pca")
DimHeatmap(wt2336blood, dims = 1, cells = 500, balanced = TRUE)

wt2336blood <- JackStraw(wt2336blood, num.replicate = 100)
wt2336blood <- ScoreJackStraw(wt2336blood, dims = 1:20)
JackStrawPlot(wt2336blood, dims = 1:15)
ElbowPlot(wt2336blood)

wt2336blood <- FindNeighbors(wt2336blood, dims = 1:15)
wt2336blood <- FindClusters(wt2336blood, resolution = 0.2)
wt2336blood <- RunTSNE(wt2336blood, dims = 1:15,seed.use = 1,perplexity=30)
DimPlot(wt2336blood, reduction = "tsne",pt.size = 2,shape.by ="condition",label = T)
FeaturePlot(wt2336blood,features = c("myh6","twist1a","vmhc"),reduction = "tsne")

wt2336blood <- RenameIdents(wt2336blood, "0" = "CM-A")
wt2336blood <- RenameIdents(wt2336blood, "5" = "CM-A")
wt2336blood <- RenameIdents(wt2336blood, "3" = "CM-V")
wt2336blood <- RenameIdents(wt2336blood, "4" = "CM-V")
wt2336blood <- RenameIdents(wt2336blood, "6" = "EP")
wt2336blood <- RenameIdents(wt2336blood, "1" = "EC")
wt2336blood <- RenameIdents(wt2336blood, "2" = "EPDC")
wt2336blood <- RenameIdents(wt2336blood, "8" = "CM-V")
wt2336blood <- RenameIdents(wt2336blood, "7" = "CM-V")
plot <- FeaturePlot(wt2336blood, features =c("myh6"),cols = c("grey", "red"),pt.size=1,reduction = "tsne")
#select.cells <- CellSelector(plot = plot)  
#reidentify the cell cluster according to the gene expression profile
select.cells <- c("B06_1dpt_MTZ_CMA.sc24", "B06_1dpt_MTZ_CMA.sc26", "B06_1dpt_MTZ_CMA.sc33",
                  "B06_1dpt_MTZ_CMA.sc58", "B06_1dpt_MTZ_CMA.sc73", "B06_1dpt_MTZ_CMA.sc74",
                  "B06_1dpt_MTZ_CMA.sc89", "B06_2dpt_CT_CMA.sc3" ,  "B06_2dpt_CT_CMA.sc5" , 
                  "B06_2dpt_CT_CMA.sc8" ,  "B06_2dpt_CT_CMA.sc10" , "B06_2dpt_CT_CMA.sc13" ,
                  "B06_2dpt_CT_CMA.sc18"  ,"B06_2dpt_CT_CMA.sc19",  "B06_2dpt_CT_CMA.sc20" ,
                  "B06_2dpt_CT_CMA.sc24" , "B06_2dpt_CT_CMA.sc30" , "B06_2dpt_CT_CMA.sc32" ,
                  "B06_2dpt_CT_CMA.sc45" , "B06_2dpt_CT_CMA.sc47" , "B06_2dpt_CT_CMV.sc55" ,
                  "B06_2dpt_MTZ_CMA.sc4" , "B06_2dpt_MTZ_CMA.sc22", "B06_2dpt_MTZ_CMA.sc49",
                  "B06_2dpt_MTZ_CMA.sc64", "B07_3dpt_CT_CMA2.sc29" ,"B07_3dpt_CT_CMA2.sc31",
                  "B06_3dpt_MTZ_CMA.sc4" , "B06_3dpt_MTZ_CMA.sc20" ,"B06_3dpt_MTZ_CMA.sc92",
                  "B06_3dpt_MTZ_CMA.sc94")
Idents(wt2336blood,cells= select.cells) <- "CM-A"
plot <- FeaturePlot(wt2336blood, features =c("twist1a"),cols = c("grey", "red"),pt.size=1,reduction = "tsne")
#select.cells2<- CellSelector(plot = plot)  
select.cells2 <- c("B01_1dpt_MTZ_V2.sc62" , "B01_1dpt_MTZ_V2.sc84",  "B01_1dpt_CT_V1.sc33"  ,
                   "B02_1dpt_CT_AV3.sc22" , "B02_1dpt_CT_AV3.sc34" , "B02_1dpt_CT_AV4.sc70" ,
                   "B02_1dpt_CT_AV4.sc79" , "B03_4dpt_CT_AV1.sc22",  "B03_4dpt_CT_AV1.sc74" ,
                   "B04_2dpt_CT_AV1.sc1"  , "B04_2dpt_CT_AV1.sc3" ,  "B04_2dpt_CT_AV2.sc34" ,
                   "B04_2dpt_MTZ_AV1.sc5" , "B04_2dpt_MTZ_AV1.sc32" ,"B04_2dpt_MTZ_AV1.sc91",
                   "B04_2dpt_MTZ_AV1.sc94","B04_2dpt_MTZ_AV1.sc95","B05_3dpt_CT_AV1.sc46",
                   "B05_3dpt_CT_AV1.sc50" , "B05_3dpt_CT_AV1.sc74" , "B05_3dpt_CT_AV2.sc49" ,
                   "B05_3dpt_MTZ_AV1.sc17" ,"B05_3dpt_MTZ_AV1.sc18" ,"B05_3dpt_MTZ_AV1.sc57",
                   "B05_3dpt_MTZ_AV1.sc63", "B05_3dpt_MTZ_AV1.sc81" ,"B05_3dpt_MTZ_AV2.sc25")
Idents(wt2336blood,cells= select.cells2) <- "EPDC"

DimPlot(wt2336blood,cells.highlight = select.cells2,reduction = "tsne")
wt2336blood <- RenameIdents(wt2336blood,"EPDC" = "eEPDC") #change order
wt2336blood <- RenameIdents(wt2336blood,"EP" = "dEP")
wt2336blood <- RenameIdents(wt2336blood,"EC" = "cEC")
wt2336blood <- RenameIdents(wt2336blood,"CM-V" = "bCM-V")
wt2336blood <- RenameIdents(wt2336blood,"CM-A" = "aCM-A")

DimPlot(wt2336blood, reduction = "tsne",pt.size = 2,shape.by ="condition",label = F,
        cols  = c("darkolivegreen4","brown3","darkgoldenrod3","dodgerblue3","darkorchid"))

# ==============================================================================================
#Fig1c, marker boxplot 
# ==============================================================================================
wt2336bloodexpress1 <- wt2336blood[["RNA"]]@counts 
wt2336bloodexpress1 <- as.data.frame(wt2336bloodexpress1)
wt2336bloodexpress1 <- t(wt2336bloodexpress1)
wt2336bloodexpress1 <- as.data.frame(wt2336bloodexpress1)
wt2336bloodexpress1 <- log2(wt2336bloodexpress1/10 + 1)
class(wt2336bloodexpress1)
dim(wt2336bloodexpress1)
table(rownames(wt2336bloodexpress1)==rownames(wt2336blood@meta.data))
Idents(wt2336blood)
table(rownames(wt2336bloodexpress1)==names(Idents(wt2336blood)))
wt2336bloodexpress1$ident <- Idents(wt2336blood)
table(wt2336bloodexpress1$ident)


library(ggplot2)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

wt2336bloodexpress1sum <- summarySE(wt2336bloodexpress1, measurevar="kdrl", groupvars="ident") #na.rm=TRUE

pdf("marker.pdf",width=1.5,height=8)
geneuse <- c("myh6","vmhc","myl7","fli1a","kdrl","tcf21","tbx18","twist1a","vim")
for(i in geneuse){
  print(i)
  wt2336bloodexpress1sum <- summarySE(wt2336bloodexpress1, measurevar=i, groupvars="ident") #na.rm=TRUE
  p <- ggplot(wt2336bloodexpress1sum,  aes_string(x="ident", y=i,fill="ident"))  +
    geom_bar(position=position_dodge(), stat="identity",width = 0.7) +
    xlab("")+ #ȥx ??ǩ
    scale_fill_manual(values = c("darkolivegreen4","brown3","darkgoldenrod3","dodgerblue3","darkorchid"))+
    #coord_fixed(ratio=1/2)+
    coord_flip()+
    scale_x_discrete(limits=rev(c("aCM-A","bCM-V","cEC","dEP","eEPDC")))+
    geom_errorbar(aes(ymin=get(i)-se, ymax=get(i)+se), # Error bars represent standard error of the mean
                  width=.1,                    # Width of the error bars
                  position=position_dodge(.9))+
    theme_bw()+ #ȥ??????ɫ
    theme(axis.text = element_blank(),#??????
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),#ȥ????????
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 25, face = "bold"),
          axis.title=element_text(size=20,face="bold"),
          legend.position = "") #ȥ??ͼ??
  
  print(p)
}
dev.off()
# =============================================================================================
#Fig1d-e, heatmap of CT/MTZ cells 
#==============================================================================================
table(wt2336blood@meta.data$condition)

wt2336bloodct <- subset(wt2336blood, condition=="CT")
table(wt2336bloodct@meta.data$condition)
wt2336bloodctmarker <- FindAllMarkers(wt2336bloodct,logfc.threshold = 0.25,test.use = "roc",only.pos = T)
library(dplyr)
wt2336bloodctmarker100 <- wt2336bloodctmarker %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
table(wt2336bloodctmarkersel$cluster)
DoHeatmap(wt2336bloodct,features = wt2336bloodctmarker100$gene,label = T,
          group.colors =  c("darkolivegreen4","brown3","darkgoldenrod3","dodgerblue3","darkorchid"))

wt2336bloodmtz <- subset(wt2336blood, condition=="MTZ")
wt2336bloodmtzmarker <- FindAllMarkers(wt2336bloodmtz,logfc.threshold = 0.25,test.use = "roc",only.pos = T)
library(dplyr)
wt2336bloodmtzmarker100 <- wt2336bloodmtzmarker %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)
DoHeatmap(wt2336bloodmtz,features = wt2336bloodmtzmarker100$gene,
          group.colors =  c("darkolivegreen4","brown3","darkgoldenrod3","dodgerblue3","darkorchid"))

save.image("~/Desktop/R code/Fig.1.RData")

