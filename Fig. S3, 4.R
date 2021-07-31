# ==============================================================================================
#Fig. S3. non-cardiomyocyte cells 
# ==============================================================================================
other <- subset(wt2336blood,idents = c("cEC","dEP","eFibro"))
dim(other@meta.data)
VlnPlot(other, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(other, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
other <- NormalizeData(other, normalization.method = "LogNormalize", scale.factor = 10000)
other <- FindVariableFeatures(other, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(other)
other <- ScaleData(other, features = all.genes)

other <- RunPCA(other, features = VariableFeatures(object = other))
ElbowPlot(other)

other <- FindNeighbors(other, dims = 1:15)
other <- FindClusters(other, resolution = 0.35)
other <- RunTSNE(other, dims = 1:15,seed.use = 5,perplexity=30)
DimPlot(other, reduction = "tsne",pt.size = 2,shape.by ="condition",label = T,
        cols = c("#F8766D","#7CAE00","#00BFC4","#C77CFF","goldenrod2"))



FeaturePlot(other, features =c("fli1a","tcf21","twist1a","mki67"),
            cols = c("grey", "red"),pt.size=1,reduction = "tsne")

FeaturePlot(other, features =c("cxcl12a"),cols = c("grey", "red"),pt.size=1,reduction = "tsne")

FeaturePlot(other, features =c("tek"),cols = c("grey", "red"),pt.size=1,reduction = "tsne")

# ==============================================================================================
##Fig. S4
# ==============================================================================================
FeaturePlot(other, features =c("slc5a1","pdgfrl","crhbp","rgs6","ugt2a1",
                               "ntd5","fstb","angpt4","gsc","pou3f1"),
            cols = c("grey", "red"),pt.size=1,reduction = "tsne",ncol = 5)

