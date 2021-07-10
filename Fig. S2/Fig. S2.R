# ==============================================================================================
#Fig. S2, a-d. analysis of cma 
# ==============================================================================================
cma <- subset(wt2336blood,idents = "aCM-A")
dim(cma@meta.data)
VlnPlot(cma, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(cma, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
cma <- NormalizeData(cma, normalization.method = "LogNormalize", scale.factor = 10000)
cma <- FindVariableFeatures(cma, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cma)
cma <- ScaleData(cma, features = all.genes) 

cma <- RunPCA(cma, features = VariableFeatures(object = cma))
ElbowPlot(cma)

cma <- FindNeighbors(cma, dims = 1:15)
cma <- FindClusters(cma, resolution = 0.25) 
DimPlot(cma, reduction = "tsne",pt.size = 2,shape.by ="condition",label = T)

cma <- RunTSNE(cma, dims = 1:15,seed.use = 1,perplexity=30)
DimPlot(cma, reduction = "tsne",pt.size = 2,shape.by ="condition",label = F,
        cols = c("sienna2","mediumorchid3","darkseagreen3","lightgoldenrod3"))
DimPlot(cma, reduction = "tsne",pt.size = 2,group.by ="timepoint",label = T)

FeaturePlot(cma, features =c("bmp4","bmp5","pcna","mki67"),
            cols = c("grey", "red"),pt.size=1,reduction = "tsne",min.cutoff = "q20")

cma <- RenameIdents(cma, "3" = "A-C4")#change order
cma <- RenameIdents(cma, "0" = "A-C3") 
cma <- RenameIdents(cma, "2" = "A-C2")
cma <- RenameIdents(cma, "1" = "A-C1")
DimPlot(cma, reduction = "tsne",pt.size = 2,shape.by ="condition",label = F,
        cols = c("mediumorchid3","darkseagreen3","sienna2","lightgoldenrod3"))

table(cma@meta.data$RNA_snn_res.0.25,cma@meta.data$condition)

ratio <- read.csv("cma condition ratioA-C1.csv",header = T)
ratio <- read.csv("cma condition ratioA-C2.csv",header = T)
ratio <- read.csv("cma condition ratioA-C3.csv",header = T)
ratio <- read.csv("cma condition ratioA-C4.csv",header = T)

ratio$per <- factor(ratio$per, levels=ratio$per)
mylabel<-paste(ratio[,2],"%") 
mylabel<-rev(mylabel)  
percent<-rev(ratio$cluster0) 
ggplot(ratio,aes(x="",y=cluster0,fill=per)) +
  geom_bar(stat = "identity",color="white") + 
  scale_fill_manual(values = c("#F8766D","#00BFC4")) +
  coord_polar(theta = "y") + 
  theme_bw()+
  theme(axis.text.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())+
  geom_text(aes(y= cumsum(percent)-percent/2, x= 1),label=mylabel)


AC1marker <- FindMarkers(cma,ident.1 = "A-C1" ,only.pos = TRUE, min.pct = 0.2, test.use = "roc") 
AC1markersel <- subset(AC1marker,power>=0.4 & avg_diff>0.6)
dim(AC1markersel)
write.csv(AC1markersel,file = "~/Desktop/R code/Fig. S2/AC1markersel.csv")

#GO boxplot
atrcmgofigure <-  read.csv("metascape_resultA-C1.csv",header = T) 
library(ggplot2)
ggplot(atrcmgofigure, aes(x=reorder(Description, -as.numeric(LogP)), y=-LogP)) + 
  geom_bar(stat="identity",fill="mediumorchid3", colour="black",width = 0.5) +
  theme_bw()+ #ȥ??????ɫ
  theme(axis.text =element_text(size=10,face="bold"),
        panel.border = element_blank(),panel.grid.major = element_blank(),#ȥ????????
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  coord_flip()

# ==============================================================================================
#Fig. S2e. cma StemID score 
# ==============================================================================================
cmaexpress <- as.data.frame(cma[["RNA"]]@counts)
dim(cmaexpress)
colnames(cmaexpress)
rownames(cmaexpress)
cmaexpress[1,1]
## RaceID2
# initialize SCseq object with transcript counts
library(RaceID)
sc <- SCseq(cmaexpress)
# filtering of expression data
sc <- filterdata(sc,mintotal = 100,minexpr = 0,minnumber = 0)
fdata <- getfdata(sc)
sc <- compdist(sc,metric="pearson")
# k-medoids clustering
sc <- clustexp(sc,clustnr=30,bootnr=50,rseed=17000,FUNcluster="kmedoids")
# compute t-SNE map
sc <- comptsne(sc,rseed=15555)
# detect outliers and redefine clusters
sc <- findoutliers(sc)
class(sc@cluster$kpart)
sc@cluster$kpart[1:10]
colnames(cmaexpress)[1:10]
as.integer(cma@active.ident[1:10])
sc@cpart <- as.integer(cma@active.ident) 
names(sc@cpart) <- names(cma@active.ident) 
table(sc@cpart)

## StemID
# initialization
ltr <- Ltree(sc)
# computation of the entropy
ltr <- compentropy(ltr)
# computation of the projections for all cells
ltr <- projcells(ltr,cthr=2,nmode=FALSE)
# computation of the projections for all cells after randomization
ltr <- projback(ltr,pdishuf=200,rseed=17000) #2000
# assembly of the lineage tree
ltr <- lineagegraph(ltr)
# determination of significant differentiation trajectories
ltr <- comppvalue(ltr)

# lineage tree showing the projections of all cells in t-SNE space
plotgraph(ltr,showCells=TRUE,scthr=0)
# lineage tree without showing the projections of all cells
plotgraph(ltr,showCells=FALSE,scthr=0)

## computing the StemID score
x <- compscore(ltr,nn=1)
x
#plotting the StemID score
plotscore(ltr,1)

x <- c('A-C1','A-C2','A-C3','A-C4') 
y <- c(0.062963317,0.054809367,0.007584146,0.00001)
df <- data.frame(x=x,y=y)
ggplot(data = df, mapping = aes(x = x, y = y,fill=x)) + geom_bar(stat= 'identity',width = 0.7)+
  scale_x_discrete(limits=c('A-C1','A-C2',"A-C3",'A-C4'))+
  scale_fill_manual(values = c("mediumorchid3","darkseagreen3","sienna2","lightgoldenrod3"))+
  theme_bw()+ 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),#ȥ????????
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size = 1),
        axis.text = element_text(size = 15),axis.title = element_blank())

# ==============================================================================================
#Fig. S2f, cma cell cycle score
# ==============================================================================================
zehuortho <- read.table("orthology.zebrafish.human.txt",sep=",",header=T)
zehuortho[1:5,]
cellcycle <- read.csv("cellcycle.csv")
cellcycle[1:5,]
table(cellcycle$cellcycle)
cellcycleg1s <- subset(cellcycle,cellcycle=="G1/S")
cellcycleg2m <- subset(cellcycle,cellcycle=="G2/M")
zeg1s <- subset(zehuortho,Human.gene.name %in% cellcycleg1s$X)
zeg2m <- subset(zehuortho,Human.gene.name %in% cellcycleg2m$X)
zecellcycle <- rbind(zeg1s,zeg2m)
zecellcycle2 <- subset(zecellcycle,zecellcycle$Gene.name %in% rownames(cma[["RNA"]]@scale.data))
zecellcycle3 <- unique(zecellcycle2$Gene.name)
table(zecellcycle3 %in% rownames(cma[["RNA"]]@scale.data))
table(zeg2m$Gene.name %in% rownames(cma[["RNA"]]@scale.data))
zeg1sgene <- intersect(zeg1s$Gene.name,zecellcycle3)
zeg2mgene <- intersect(zeg2m$Gene.name,zecellcycle3)

expressg1s <- cma[["RNA"]]@counts[zeg1sgene,]
expressg1s <- as.matrix(expressg1s)
dim(expressg1s)
colMeans(expressg1s,na.rm = T)
table(colMeans(expressg1s,na.rm = T) > 2)

expressg2m <- cma[["RNA"]]@counts[zeg2mgene,]
expressg2m <- as.matrix(expressg2m)
dim(expressg2m)
colMeans(expressg2m,na.rm = T)
table(colMeans(expressg2m,na.rm = T) > 2)

cyclescore <- data.frame(colMeans(expressg1s,na.rm = T),colMeans(expressg2m,na.rm = T))
dim(cyclescore)
colnames(cyclescore) <- c("g1s","g2m")
prolicell <- rownames(subset(cyclescore,g1s>1.5  | g2m>1.5)) 
length(prolicell)

# by ident
scoreA-C1 <- length(intersect(prolicell,WhichCells(cma,ident = "A-C1"))) / length(WhichCells(cma,ident = "A-C1"))
scoreA-C2 <- length(intersect(prolicell,WhichCells(cma,ident = "A-C2"))) / length(WhichCells(cma,ident = "A-C2"))
scoreA-C3 <- length(intersect(prolicell,WhichCells(cma,ident = "A-C3"))) / length(WhichCells(cma,ident = "A-C3"))
scoreA-C4 <- length(intersect(prolicell,WhichCells(cma,ident = "A-C4"))) / length(WhichCells(cma,ident = "A-C4"))


x <- c('A-C1','A-C2','A-C3','A-C4') 
y <- c(scoreA-C1,scoreA-C2,scoreA-C3,scoreA-C4)
df <- data.frame(x=x,y=y)
ggplot(data = df, mapping = aes(x = x, y = y,fill=x)) + geom_bar(stat= 'identity',width = 0.7)+
  scale_x_discrete(limits=c('A-C1','A-C2','A-C3','A-C4'))+
  scale_fill_manual(values = c("mediumorchid3","darkseagreen3","sienna2","lightgoldenrod3"))+
  theme_bw()+ 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size = 1),
        axis.text = element_text(size = 15),axis.title = element_blank())

# ==============================================================================================
#Fig. S2g, cma transdifferentiation pseudotime 
# ==============================================================================================
cmatrans <- subset(wt2336blood,cells=intersect(WhichCells(cma,idents = "A-C1"),WhichCells(cma,expression = condition=="MTZ")))
dim(cmatrans@meta.data)
cmatrans@meta.data$av <- rep("amtz",117)

cmvtrans <- subset(cmv,cells = WhichCells(cmv,expression = condition=="MTZ"))
table(cmvtrans@meta.data$RNA_snn_res.0.3)
cmvtrans@meta.data$av <- cmvtrans@meta.data$RNA_snn_res.0.3

cmtrans <- merge(cmvtrans,cmatrans)
dim(cmtrans@meta.data)
table(cmtrans@meta.data$av)

library(monocle)
dim(cmtrans[["RNA"]]@counts)
class(cmtrans[["RNA"]]@counts)
cmtransexpress2 <- as.matrix(cmtrans[["RNA"]]@counts)
dim(cmtransexpress2)
dim(cmtrans@meta.data)
cmtrans@meta.data[1:6,]
cmtransexpress2[1:6,1:4]
cmtranspd <- cmtrans@meta.data
colnames(cmtranspd)
cmtranspd$modi <- as.character(cmtranspd$timepoint) 
cmtranspd[rownames(subset(cmtrans@meta.data,condition=="Con")),"modi"]= rep("uninjure",140)
dim(cmtranspd)
table(cmtranspd$av,cmtranspd$condition)

rownames(cmtranspd) <- rownames(cmtrans@meta.data)
table(rownames(cmtranspd)==colnames(cmtransexpress2))
cmtrans41 <- rownames(cmtransexpress2)
cmtransfd <- as.data.frame(cmtrans41)
rownames(cmtransfd) <- cmtransfd[,1]
colnames(cmtransfd) <- "gene_short_name"
table(rownames(cmtransfd)==rownames(cmtransexpress2))
pd <- new("AnnotatedDataFrame",data=cmtranspd)
fd <- new("AnnotatedDataFrame",data=cmtransfd)
pseucmtrans <- newCellDataSet(as.matrix(cmtransexpress2),phenoData = pd,featureData = fd,expressionFamily = negbinomial())
pseucmtrans <- estimateSizeFactors(pseucmtrans)
pseucmtrans <- estimateDispersions(pseucmtrans)

pseucmtrans <- detectGenes(pseucmtrans,min_expr = 0.1)
pseucmtransexpressed_genes <- row.names(cmtrans[["RNA"]]@scale.data)
length(pseucmtransexpressed_genes)

pseucmtransdiff_test_res <- differentialGeneTest(pseucmtrans[pseucmtransexpressed_genes,],cores = 4,
                                                  fullModelFormulaStr = "~av")

pseucmtransordering_genes <- row.names (subset(pseucmtransdiff_test_res, 
                                               qval < 0.0000001)) 
length(pseucmtransordering_genes)
pseucmtrans <- setOrderingFilter(pseucmtrans, pseucmtransordering_genes)


pseucmtrans <- reduceDimension(pseucmtrans, max_components = 2,reduction_method = 'DDRTree',norm_method = "none")
pseucmtrans <- orderCells(pseucmtrans)

plot_cell_trajectory(pseucmtrans, color_by = "av")+
  scale_color_manual(values =c("#FFC125","#264BCC","#7FA520", "#0ABECC","#CD3333","mediumorchid3")) 

plot_cell_trajectory(pseucmtrans, color_by = "RNA_snn_res.0.3")
plot_cell_trajectory(pseucmtrans, color_by = "timepoint")
plot_cell_trajectory(pseucmtrans, color_by = "Pseudotime")
plot_cell_trajectory(pseucmtrans, color_by = "State")



