library(dplyr)
library(Seurat)
library(patchwork)
rm(list = ls())
# Load the dataset
mtx_dir_ad1 = "D:/MS/Job/MGH/mtx/AD1_AD2/"
mtx_dir_ad3 = "D:/MS/Job/MGH/mtx/AD3_AD4/"

ad1.data <- Read10X(data.dir = mtx_dir_ad1)
ad1 <- CreateSeuratObject(counts = ad1.data, project = "ad1", min.cells = 3, min.features = 200)
ad1 <- RenameCells(ad1, add.cell.id = "AD1_AD2")
ad1[["percent.mt"]] <- PercentageFeatureSet(ad1, pattern = "^MT-")

ad3.data <- Read10X(data.dir = mtx_dir_ad3)
ad3 <- CreateSeuratObject(counts = ad3.data, project = "ad3", min.cells = 3, min.features = 200)
ad3 <- RenameCells(ad3, add.cell.id = "AD3_AD4")
ad3[["percent.mt"]] <- PercentageFeatureSet(ad3, pattern = "^MT-")

ad_13 <- merge(ad1, ad3)
table(ad_13$orig.ident)

ad_13 <- subset(ad_13, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 7)
ad_13.list <- SplitObject(ad_13, split.by="orig.ident")


ad_13.list <- lapply(X = ad_13.list, FUN = function(x) {
  # Normalization
  x <- NormalizeData(x, scale.factor = 10000)
  # Feature selection
  x <- FindVariableFeatures(x, x.low.cutoff=0.0125, x.high.cutoff=3,y.cutoff=0.5)
  # Scaling the data
  x <- ScaleData(x)
})

features <- SelectIntegrationFeatures(object.list = ad_13.list)
ad_13.list <- lapply(X = ad_13.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

immune.anchors <-  FindIntegrationAnchors(object.list = ad_13.list, anchor.features = features, reduction = "rpca")
immune.combined <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)

immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)























