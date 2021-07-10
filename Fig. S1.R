# ==============================================================================================
#Fig. S1 a-d
# ==============================================================================================
load("~/Desktop/R code/Fig.1.RData")

DimPlot(wt2336blood, reduction = "tsne",pt.size = 1,cells = WhichCells(wt2336bloodct,expression = timepoint=="1dpt"),
        cols  = c("darkolivegreen4","brown3","darkgoldenrod3","dodgerblue3","darkorchid"))
DimPlot(wt2336blood, reduction = "tsne",pt.size = 1,cells = WhichCells(wt2336bloodct,expression = timepoint=="2dpt"),
        cols  = c("darkolivegreen4","brown3","darkgoldenrod3","dodgerblue3","darkorchid"))
DimPlot(wt2336blood, reduction = "tsne",pt.size = 1,cells = WhichCells(wt2336bloodct,expression = timepoint=="3dpt"),
        cols  = c("darkolivegreen4","brown3","darkgoldenrod3","dodgerblue3","darkorchid"))
DimPlot(wt2336blood, reduction = "tsne",pt.size = 1,cells = WhichCells(wt2336bloodct,expression = timepoint=="4dpt"),
        cols  = c("darkolivegreen4","brown3","darkgoldenrod3","dodgerblue3","darkorchid"))

DimPlot(wt2336blood, reduction = "tsne",pt.size = 1,cells = WhichCells(wt2336bloodmtz,expression = timepoint=="1dpt"),
        cols  = c("darkolivegreen4","brown3","darkgoldenrod3","dodgerblue3","darkorchid"))
DimPlot(wt2336blood, reduction = "tsne",pt.size = 1,cells = WhichCells(wt2336bloodmtz,expression = timepoint=="2dpt"),
        cols  = c("darkolivegreen4","brown3","darkgoldenrod3","dodgerblue3","darkorchid"))
DimPlot(wt2336blood, reduction = "tsne",pt.size = 1,cells = WhichCells(wt2336bloodmtz,expression = timepoint=="3dpt"),
        cols  = c("darkolivegreen4","brown3","darkgoldenrod3","dodgerblue3","darkorchid"))
DimPlot(wt2336blood, reduction = "tsne",pt.size = 1,cells = WhichCells(wt2336bloodmtz,expression = timepoint=="4dpt"),
        cols  = c("darkolivegreen4","brown3","darkgoldenrod3","dodgerblue3","darkorchid"))


hist(as.numeric(wt2336blood@meta.data$nFeature_RNA/1000),breaks = 100,freq = T,col = "#74A6DA",
     border = NA,xlim = c(1,5))
hist(log10(as.numeric(wt2336blood@meta.data$nCount_RNA)),breaks = 100,freq = T,col = "#74A6DA",
     border = NA,xlim = c(3.5,5.5))
# ==============================================================================================
#Fig. S1 e
# ==============================================================================================
dim(wt2336blood@meta.data)
table(wt2336blood@active.ident)
expressallnewlog <- log2((expresswt2336[,rownames(wt2336blood@meta.data)]/10)+1) 
dim(expressallnewlog)
expressallnewlog[1:5,1:5]

wt2336bloodcluster <- as.data.frame(t(expressallnewlog))
wt2336bloodcluster$ident <- as.character(wt2336blood@active.ident) 
dim(wt2336bloodcluster)
epicell <- WhichCells(wt2336blood,ident = "dEP")
endocell <- WhichCells(wt2336blood,ident = "cEC")
fibcell <- WhichCells(wt2336blood,ident="eFibro")
cmacell <- WhichCells(wt2336blood,ident="aCM-A")
cmvcell <- WhichCells(wt2336blood,ident="bCM-V")
length(epicell)
length(endocell)
length(fibcell)
length(cmacell)
length(cmvcell)

table(c(cmacell,cmvcell,epicell,endocell,fibcell) %in% colnames(expressallnewlog))
wt2336bloodcluster <- wt2336bloodcluster[c(cmacell,cmvcell,endocell,epicell,fibcell),]  
wt2336bloodcluster$num <- as.factor(c(1:1581))
dim(wt2336bloodcluster)
colnames(wt2336bloodcluster)
rownames(wt2336bloodcluster)

colors = c(rep("darkolivegreen4",length(cmacell)),rep("brown3",length(cmvcell)),
           rep("darkgoldenrod3",length(endocell)),rep("dodgerblue3",length(epicell)), rep("darkorchid",length(fibcell)))
length(colors)

ggplot(wt2336bloodcluster,aes(num,actb2))+geom_histogram(stat = 'identity',colour=colors)+ #6*5  #aes vs aes_string
  theme(axis.line = element_line(size=0.8, colour = "black"),axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) 
