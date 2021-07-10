library(Seurat)
library(ggplot2)
library(monocle)
setwd("~/Desktop/R code/")
load("~/Desktop/R code/Fig.2.RData")
# =============================================================================================
#ang mutant seq, cma and cmv
#==============================================================================================
expressang576 <- read.csv("angexpress576.csv",stringsAsFactors=F)
dim(expressang576)
rownames(expressang576) <- expressang576[,1]
rownames(expressang576)
grep("ERCC",rownames(expressang576),value = T)
expressang576 <- expressang576[setdiff(rownames(expressang576),grep("ERCC",rownames(expressang576),value = T)),] #ȥ??ERCC
grep("mt-",rownames(expressang576),value = T)
colnames(expressang576)
expressang576[1,1]
expressang576 <- expressang576[,-1]
table(rownames(expresswt2336)==rownames(expressang576))

metadataang576 <- read.csv("angmeta576.csv")
rownames(metadataang576) <- metadataang576[,1]
colnames(metadataang576)
table(rownames(metadataang576)==colnames(expressang576))

ang576 <- CreateSeuratObject(counts = expressang576, project = "ang576", min.cells = 3,min.features = 1000) 
ang576 <- AddMetaData(ang576,metadata = metadataang576 )
dim(ang576@meta.data)
table(ang576@meta.data$VA)

# =============================================================================================
#Fig 5a, cmvall, combine angcmv and wtcmv for downstream analysis
#==============================================================================================
angcmv <- subset(ang576, subset = VA == "CMV")
cmvall <- merge(angcmv,cmv)
dim(cmvall@meta.data)
table(cmvall@meta.data$condition,cmvall@meta.data$run)
cmvall@meta.data$modi <- as.character(cmvall@meta.data$timepoint) 
cmvall@meta.data[WhichCells(object = cmvall, expression = condition =="CT"),"modi"]= rep("uninjure",212)
cmvall@meta.data$modi2 <- cmvall@meta.data$modi 
cmvall@meta.data[WhichCells(object = cmvall, expression = run=="B08"),"modi2"]= rep("ang",244)
cmvall@meta.data[WhichCells(object = cmvall, expression = run !="B08" & condition=="MTZ"),"modi2"]= rep("WTMTZ",298)

dim(cmvall@meta.data)
table(cmvall@meta.data$modi2)
cmvall[["percent.mt"]] <- PercentageFeatureSet(cmvall, pattern = "^mt-")
VlnPlot(cmvall, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(cmvall, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
cmvall <- subset(cmvall, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & nCount_RNA < 200000)

cmvall <- NormalizeData(cmvall, normalization.method = "LogNormalize", scale.factor = 10000)
cmvall <- FindVariableFeatures(cmvall, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(cmvall), 10)
plot1 <- VariableFeaturePlot(cmvall)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

all.genes <- rownames(cmvall)
cmvall <- ScaleData(cmvall, features = all.genes) 
cmvall <- RunPCA(cmvall, features = VariableFeatures(object = cmvall))
print(cmvall[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(cmvall, dims = 1:2, reduction = "pca")
DimPlot(cmvall, reduction = "pca")
DimHeatmap(cmvall, dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(cmvall)

cmvall <- FindNeighbors(cmvall, dims = 1:15)
cmvall <- FindClusters(cmvall, resolution = 0.5)
cmvall <- RunTSNE(cmvall, dims = 1:15,seed.use = 2,perplexity=25)
DimPlot(cmvall, reduction = "tsne",pt.size = 2,shape.by ="modi2",label = T)

FeaturePlot(cmvall,reduction = "tsne",features = c("fzd2","rac1a","ctsk","ctsz"),
            cols = c("grey","red"),min.cutoff = c("q10","q10","q50","q50"))

# =============================================================================================
#Fig. 5b-d cmvallmtz pseudotime
#==============================================================================================
cmvallmtz <- subset(cmvall,subset=condition=="MTZ")
dim(cmvallmtz@meta.data)
table(cmvallmtz@meta.data$RNA_snn_res.0.5)
table(cmvall@meta.data$RNA_snn_res.0.5,cmvall@meta.data$condition)

library(monocle)
dim(cmvallmtz[["RNA"]]@counts)
class(cmvallmtz[["RNA"]]@counts)
cmvallmtzexpress2 <- as.matrix(cmvallmtz[["RNA"]]@counts)
dim(cmvallmtzexpress2)
dim(cmvallmtz@meta.data)
cmvallmtz@meta.data[1:6,]
cmvallmtzexpress2[1:6,1:4]
cmvallmtzpd <- cmvallmtz@meta.data
colnames(cmvallmtzpd)
cmvallmtzpd$modi <- as.character(cmvallmtzpd$timepoint) 
cmvallmtzpd[rownames(subset(cmvallmtz@meta.data,condition=="Con")),"modi"]= rep("uninjure",140)
dim(cmvallmtzpd)
table(cmvallmtzpd$modi2,cmvallmtzpd$condition)
table(cmvallmtzpd$RNA_snn_res.0.5)

rownames(cmvallmtzpd) <- rownames(cmvallmtz@meta.data)
table(rownames(cmvallmtzpd)==colnames(cmvallmtzexpress2))
cmvallmtz41 <- rownames(cmvallmtzexpress2)
cmvallmtzfd <- as.data.frame(cmvallmtz41)
rownames(cmvallmtzfd) <- cmvallmtzfd[,1]
colnames(cmvallmtzfd) <- "gene_short_name"
table(rownames(cmvallmtzfd)==rownames(cmvallmtzexpress2))
pd <- new("AnnotatedDataFrame",data=cmvallmtzpd)
fd <- new("AnnotatedDataFrame",data=cmvallmtzfd)
pseucmvallmtz <- newCellDataSet(as.matrix(cmvallmtzexpress2),phenoData = pd,featureData = fd,expressionFamily = negbinomial())
pseucmvallmtz <- estimateSizeFactors(pseucmvallmtz)
pseucmvallmtz <- estimateDispersions(pseucmvallmtz)

pseucmvallmtz <- detectGenes(pseucmvallmtz,min_expr = 0.1)
pseucmvallmtzexpressed_genes <- row.names(cmvallmtz[["RNA"]]@scale.data)
length(pseucmvallmtzexpressed_genes)

pseucmvallmtzdiff_test_res <- differentialGeneTest(pseucmvallmtz[pseucmvallmtzexpressed_genes,],cores = 4,
                                                    fullModelFormulaStr = "~RNA_snn_res.0.5")

pseucmvallmtzordering_genes <- row.names (subset(pseucmvallmtzdiff_test_res, qval < 0.00000000000000001)) 
length(pseucmvallmtzordering_genes)
pseucmvallmtz <- setOrderingFilter(pseucmvallmtz, pseucmvallmtzordering_genes)

pseucmvallmtz <- reduceDimension(pseucmvallmtz, max_components = 10,reduction_method = 'DDRTree',norm_method = "none",cores = 4)
pseucmvallmtz <- orderCells(pseucmvallmtz)

plot_cell_trajectory(pseucmvallmtz, color_by = "Pseudotime", show_branch_points = F)
plot_cell_trajectory(pseucmvallmtz, color_by = "RNA_snn_res.0.5", show_branch_points = F)
plot_cell_trajectory(pseucmvallmtz, color_by = "modi2", show_branch_points = T) + 
  scale_color_manual(values= c("#87479D","#F38401")) 

table(cmvallmtz@meta.data$modi2)
angmtzcell <- WhichCells(object = cmvallmtz, expression = modi2== "ang")
wtmtzcell <- WhichCells(object = cmvallmtz, expression = modi2== "WTMTZ")

monocle::plot_genes_in_pseudotime(pseucmvallmtz["vmhc",],color_by = "modi2",cell_size = 1.5)+scale_color_manual(values= c("#87479D","#F38401"))+
  scale_y_log10()+
  theme(axis.line.x  = element_line(size = 0.6,colour = "black"),axis.line.y  = element_line(size = 0.6,colour = "black"),
        axis.text.x = element_text(size = 15,colour = "black"),axis.text.y = element_text(size = 15,colour = "black"))#+facet_wrap(~modi2, nrow = 1)
monocle::plot_genes_in_pseudotime(pseucmvallmtz["vmhc",angmtzcell],color_by = "modi2",cell_size = 1.5)+ scale_color_manual(values= c("#87479D"))+
  scale_y_log10()+
  theme(axis.line.x  = element_line(size = 0.6,colour = "black"),axis.line.y  = element_line(size = 0.6,colour = "black"),
        axis.text.x = element_text(size = 15,colour = "black"),axis.text.y = element_text(size = 15,colour = "black"))
monocle::plot_genes_in_pseudotime(pseucmvallmtz["vmhc",wtmtzcell],color_by = "modi2",cell_size = 1.5) + scale_color_manual(values= c("#F38401"))+
  theme(axis.line.x  = element_line(size = 0.6,colour = "black"),axis.line.y  = element_line(size = 0.6,colour = "black"),
        axis.text.x = element_text(size = 15,colour = "black"),axis.text.y = element_text(size = 15,colour = "black"))

monocle::plot_genes_in_pseudotime(pseucmvallmtz["tnnt2a",],color_by = "modi2",cell_size = 1.5)+scale_color_manual(values= c("#87479D","#F38401"))+
  scale_y_log10()+
  theme(axis.line.x  = element_line(size = 0.6,colour = "black"),axis.line.y  = element_line(size = 0.6,colour = "black"),
        axis.text.x = element_text(size = 15,colour = "black"),axis.text.y = element_text(size = 15,colour = "black"))#+facet_wrap(~modi2, nrow = 1)
monocle::plot_genes_in_pseudotime(pseucmvallmtz["tnnt2a",angmtzcell],color_by = "modi2",cell_size = 1.5)+ scale_color_manual(values= c("#87479D"))+
  scale_y_log10()+
  theme(axis.line.x  = element_line(size = 0.6,colour = "black"),axis.line.y  = element_line(size = 0.6,colour = "black"),
        axis.text.x = element_text(size = 15,colour = "black"),axis.text.y = element_text(size = 15,colour = "black"))
monocle::plot_genes_in_pseudotime(pseucmvallmtz["tnnt2a",wtmtzcell],color_by = "modi2",cell_size = 1.5) + scale_color_manual(values= c("#F38401"))+
  theme(axis.line.x  = element_line(size = 0.6,colour = "black"),axis.line.y  = element_line(size = 0.6,colour = "black"),
        axis.text.x = element_text(size = 15,colour = "black"),axis.text.y = element_text(size = 15,colour = "black"))


monocle::plot_genes_in_pseudotime(pseucmvallmtz["tnnt2b",],color_by = "modi2",cell_size = 1.5)+scale_color_manual(values= c("#87479D","#F38401"))+
  scale_y_log10()+
  theme(axis.line.x  = element_line(size = 0.6,colour = "black"),axis.line.y  = element_line(size = 0.6,colour = "black"),
        axis.text.x = element_text(size = 15,colour = "black"),axis.text.y = element_text(size = 15,colour = "black"))#+facet_wrap(~modi2, nrow = 1)
monocle::plot_genes_in_pseudotime(pseucmvallmtz["tnnt2b",angmtzcell],color_by = "modi2",cell_size = 1.5)+ scale_color_manual(values= c("#87479D"))+
  scale_y_log10()+
  theme(axis.line.x  = element_line(size = 0.6,colour = "black"),axis.line.y  = element_line(size = 0.6,colour = "black"),
        axis.text.x = element_text(size = 15,colour = "black"),axis.text.y = element_text(size = 15,colour = "black"))
monocle::plot_genes_in_pseudotime(pseucmvallmtz["tnnt2b",wtmtzcell],color_by = "modi2",cell_size = 1.5) + scale_color_manual(values= c("#F38401"))+
  theme(axis.line.x  = element_line(size = 0.6,colour = "black"),axis.line.y  = element_line(size = 0.6,colour = "black"),
        axis.text.x = element_text(size = 15,colour = "black"),axis.text.y = element_text(size = 15,colour = "black"))





