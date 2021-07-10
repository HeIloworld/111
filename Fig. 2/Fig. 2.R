library(Seurat)
library(ggplot2)
library(monocle)
setwd("~/Desktop/R code/")
# ==============================================================================================
#Fig 2 a-b 
# ==============================================================================================
load("~/Desktop/R code/Fig.1.RData")
setwd("~/Desktop/R code/Fig. 2/")

cmv <- subset(wt2336blood,idents = "bCM-V")
dim(cmv@meta.data)
VlnPlot(cmv, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(cmv, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
cmv <- NormalizeData(cmv, normalization.method = "LogNormalize", scale.factor = 10000)
cmv <- FindVariableFeatures(cmv, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cmv)
cmv <- ScaleData(cmv, features = all.genes) #vars.to.regress = "percent.mt"

cmv <- RunPCA(cmv, features = VariableFeatures(object = cmv))
ElbowPlot(cmv)

cmv <- FindNeighbors(cmv, dims = 1:15)
cmv <- FindClusters(cmv, resolution = 0.3) 
DimPlot(cmv, reduction = "tsne",pt.size = 2,shape.by ="condition",label = T)

cmv <- RenameIdents(cmv, "2" = "V-C5")
cmv <- RenameIdents(cmv, "4" = "V-C4")
cmv <- RenameIdents(cmv, "0" = "V-C3")
cmv <- RenameIdents(cmv, "3" = "V-C2")
cmv <- RenameIdents(cmv, "1" = "V-C1")

cmv <- RunTSNE(cmv, dims = 1:15,seed.use = 2,perplexity=20)
DimPlot(cmv, reduction = "tsne",pt.size = 2,shape.by ="condition",label = T,
        cols = c("#264BCC","#0ABECC","#FFC125","#CD3333","#7FA520"))

FeaturePlot(cmv, features =c("myl4","tnnt2b","hand2","nkx2.5","cxcr4b","twist1a"),cols = c("grey", "red"),
            pt.size=1,reduction = "tsne")


VC1marker <- FindMarkers(cmv,ident.1 = "V-C1" ,only.pos = TRUE, min.pct = 0.2, test.use = "roc") 
VC1markersel <- subset(VC1marker,power>=0.4 & avg_diff>0.6)

VC2marker <- FindMarkers(cmv,ident.1 = "V-C2" ,only.pos = TRUE, min.pct = 0.2, test.use = "roc") 
VC2markersel <- subset(VC2marker,power>=0.4 & avg_diff>0.6)

VC3marker <- FindMarkers(cmv,ident.1 = "V-C3" ,only.pos = TRUE, min.pct = 0.2, test.use = "roc") 
VC3markersel <- subset(VC3marker,power>=0.4 & avg_diff>0.6)

VC4marker <- FindMarkers(cmv,ident.1 = "V-C4" ,only.pos = TRUE, min.pct = 0.2, test.use = "roc") 
VC4markersel <- subset(VC4marker,power>=0.4 & avg_diff>0.6)

VC5marker <- FindMarkers(cmv,ident.1 = "V-C5" ,only.pos = TRUE, min.pct = 0.2, test.use = "roc") 
VC5markersel <- subset(VC5marker,power>=0.4 & avg_diff>0.6)

# ==============================================================================================
#Fig 2c, cmv cell cycle score
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
zecellcycle2 <- subset(zecellcycle,zecellcycle$Gene.name %in% rownames(cmv[["RNA"]]@scale.data))
zecellcycle3 <- unique(zecellcycle2$Gene.name)
table(zecellcycle3 %in% rownames(cmv[["RNA"]]@scale.data))
table(zeg2m$Gene.name %in% rownames(cmv[["RNA"]]@scale.data))
zeg1sgene <- intersect(zeg1s$Gene.name,zecellcycle3)
zeg2mgene <- intersect(zeg2m$Gene.name,zecellcycle3)

#cell cycle score
expressg1s <- cmv[["RNA"]]@counts[zeg1sgene,]
expressg1s <- as.matrix(expressg1s)
dim(expressg1s)
colMeans(expressg1s,na.rm = T)
table(colMeans(expressg1s,na.rm = T) > 2)

expressg2m <- cmv[["RNA"]]@counts[zeg2mgene,]
expressg2m <- as.matrix(expressg2m)
dim(expressg2m)
colMeans(expressg2m,na.rm = T)
table(colMeans(expressg2m,na.rm = T) > 2)

cyclescore <- data.frame(colMeans(expressg1s,na.rm = T),colMeans(expressg2m,na.rm = T))
dim(cyclescore)
colnames(cyclescore) <- c("g1s","g2m")
rownames(cyclescore)
prolicell <- rownames(subset(cyclescore,g1s>1.5  | g2m>1.5)) #2 1.5
length(prolicell)

# score by cluster
scoreVC1 <- length(intersect(prolicell,WhichCells(cmv,ident = "V-C1"))) / length(WhichCells(cmv,ident = "V-C1"))
scoreVC2 <- length(intersect(prolicell,WhichCells(cmv,ident = "V-C2"))) / length(WhichCells(cmv,ident = "V-C2"))
scoreVC3 <- length(intersect(prolicell,WhichCells(cmv,ident = "V-C3"))) / length(WhichCells(cmv,ident = "V-C3"))
scoreVC4 <- length(intersect(prolicell,WhichCells(cmv,ident = "V-C4"))) / length(WhichCells(cmv,ident = "V-C4"))
scoreVC5 <- length(intersect(prolicell,WhichCells(cmv,ident = "V-C5"))) / length(WhichCells(cmv,ident = "V-C5"))

x <- c('scoreVC1','scoreVC2','scoreVC3','scoreVC4','scoreVC5') 
y <- c(scoreVC1,scoreVC2,scoreVC3,scoreVC4,scoreVC5)
df <- data.frame(x=x,y=y)
ggplot(data = df, mapping = aes(x = x, y = y,fill=x)) + geom_bar(stat= 'identity',width = 0.7)+
  scale_x_discrete(limits=c('scoreVC1','scoreVC2','scoreVC3','scoreVC4','scoreVC5'))+
  scale_fill_manual(values = c("#264BCC","#0ABECC","#FFC125","#CD3333","#7FA520"))+
  theme_bw()+ 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),#ȥ????????
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size = 1),
        axis.text = element_text(size = 15),axis.title = element_blank())

# ==============================================================================================
#Fig 2d-e cmv-mtz pseudotime
# ==============================================================================================
cmvmtz <- subset(cmv,condition=="MTZ")
library(monocle)
dim(cmvmtz[["RNA"]]@counts)
class(cmvmtz[["RNA"]]@counts)
cmvmtzexpress2 <- as.matrix(cmvmtz[["RNA"]]@counts)
dim(cmvmtzexpress2)
dim(cmvmtz@meta.data)
cmvmtz@meta.data[1:6,]
cmvmtzexpress2[1:6,1:4]
cmvmtzpd <- cmvmtz@meta.data
cmvmtzpd$ident <- as.character(cmvmtz@active.ident)
colnames(cmvmtzpd)
dim(cmvmtzpd)

rownames(cmvmtzpd) <- rownames(cmvmtz@meta.data)
table(rownames(cmvmtzpd)==colnames(cmvmtzexpress2))
cmvmtz41 <- rownames(cmvmtzexpress2)
cmvmtzfd <- as.data.frame(cmvmtz41)
rownames(cmvmtzfd) <- cmvmtzfd[,1]
colnames(cmvmtzfd) <- "gene_short_name"
table(rownames(cmvmtzfd)==rownames(cmvmtzexpress2))
pd <- new("AnnotatedDataFrame",data=cmvmtzpd)
fd <- new("AnnotatedDataFrame",data=cmvmtzfd)
pseucmvmtz <- newCellDataSet(as.matrix(cmvmtzexpress2),phenoData = pd,featureData = fd,expressionFamily = negbinomial())
pseucmvmtz <- estimateSizeFactors(pseucmvmtz)
pseucmvmtz <- estimateDispersions(pseucmvmtz)

pseucmvmtz <- detectGenes(pseucmvmtz,min_expr = 0.1)
pseucmvmtzexpressed_genes <- row.names(cmvmtz[["RNA"]]@scale.data)
length(pseucmvmtzexpressed_genes)
pseucmvmtzdiff_test_res <- differentialGeneTest(pseucmvmtz[pseucmvmtzexpressed_genes,],cores = 4,
                                                 fullModelFormulaStr = "~ident")

pseucmvmtzordering_genes <- row.names (subset(pseucmvmtzdiff_test_res, 
                                              qval < 0.000000000000001)) 
length(pseucmvmtzordering_genes)
pseucmvmtz <- setOrderingFilter(pseucmvmtz, pseucmvmtzordering_genes)

pseucmvmtz <- reduceDimension(pseucmvmtz, max_components = 2,reduction_method = 'DDRTree',norm_method = "none")
pseucmvmtz <- orderCells(pseucmvmtz)

plot_cell_trajectory(pseucmvmtz, color_by = "ident",cell_size = 2) + 
  scale_color_manual(values =c("#264BCC","#0ABECC","#FFC125","#CD3333","#7FA520")) 

plot_cell_trajectory(pseucmvmtz, color_by = "State",cell_size = 2)
pseucmvmtz <- orderCells(pseucmvmtz,root_state = 3) 
plot_cell_trajectory(pseucmvmtz, color_by = "Pseudotime",cell_size = 2)

plot_pseudotime_heatmap(pseucmvmtz[intersect(unique(c(rownames(VC1markersel)[1:50],rownames(VC2markersel)[1:50],
                                                      rownames(VC3markersel)[1:50],rownames(VC4markersel)[1:50],
                                                      rownames(VC5markersel)[1:50],
                                                      "sgcd","qki2","ctsc","cxcr4b","csf3b","fbln5","lpar1","dact2",
                                                      "klf2a","tbx5a","hand2")), pseucmvmtzexpressed_genes),],
                        cores = 4,show_rownames = F,cluster_rows = T,num_clusters = 4)

# ==============================================================================================
#go mtz pseudotime, different stage
# ==============================================================================================
cmvmtzlog <- log2(cmvmtz[["RNA"]]@counts/10 + 1)
class(cmvmtzlog)
dim(cmvmtzlog)

tissueregeneration <- read.table("tissue regeneration.txt",sep = "\t",header = F)
tissueregenerationgene <- unique(tissueregeneration$V15) 
geneuse <- intersect(tissueregenerationgene,rownames(cmvmtz[["RNA"]]@counts))

mesodermdevelopment <- read.table("mesoderm development.txt",sep = "\t",header = F)
mesodermdevelopmentgene <- unique(mesodermdevelopment$V15) 
geneuse <- intersect(mesodermdevelopmentgene,rownames(cmvmtz[["RNA"]]@counts))

cellmigration <- read.table("cell migration.txt",sep = "\t",header = F)
cellmigrationgene <- unique(cellmigration$V15) 
geneuse <- intersect(cellmigrationgene,rownames(cmvmtz[["RNA"]]@counts))

ROSmetabolicprocess <- read.table("reactive oxygen species metabolic process.txt",sep = "\t",header = F)
ROSmetabolicprocessgene <- unique(ROSmetabolicprocess$V15) 
geneuse <- intersect(ROSmetabolicprocessgene,rownames(cmvmtz[["RNA"]]@counts))

proteolysis <- read.table("proteolysis.txt",sep = "\t",header = F)
proteolysisgene <- unique(proteolysis$V15) 
geneuse <- intersect(proteolysisgene,rownames(cmvmtz[["RNA"]]@counts))

cmvmtzexpress2go1 <- cmvmtz[["RNA"]]@scale.data[geneuse,]
dim(cmvmtzexpress2go1)
cmvmtzexpress2go1sum <- colSums(cmvmtzexpress2go1) #try colsums or colMeans
table(names(cmvmtzexpress2go1sum)==rownames(pseucmvmtz@phenoData@data))
go1fin <- data.frame(cmvmtzexpress2go1sum,pseucmvmtz@phenoData@data)
go1fin$Pseudotimefac <- as.numeric(factor(go1fin$Pseudotime))
dim(go1fin)
ggplot(go1fin, aes(x=Pseudotimefac, y=cmvmtzexpress2go1sum))  + geom_point(aes(colour = Pseudotime) )+ 
  theme_bw()+ 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black",size = 1),
        axis.text = element_text(size = 15),axis.title = element_blank())+ 
  stat_smooth(method = loess,se=T,color="grey15",size=1.5)

save.image("~/Desktop/R code/Fig.2.RData")
