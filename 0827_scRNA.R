library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
rm(list = ls())
assays <- dir("./mtx/")
dir <- paste0("./mtx/", assays)
samples_name = c("AD1_AD2","AD3_AD4","AD5_AD6","Ct1_Ct2","Ct3_Ct4","Ct5_Ct6")

scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = samples_name[i],
                                       min_cells = 3, min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])
  if(T){
    scRNAlist[[i]][["percent.mt"]]<- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  }
}

names(scRNAlist) <- samples_name
scRNA <- merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])

plot.features =c("nFeature_RNA","nCount_RNA","percent.mt")
plots = list()
for(i in seq_along(plot.features)){
  plots[[i]] = VlnPlot(scRNA, group.by = "orig.ident", pt.size = 0,
                       features = plot.features[i]) + NoLegend()
}
violin <- wrap_plots(plots = plots, nrow = 2)
dir.create("QC")
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 9, height = 8)

scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 7)
plots = list()
for(i in seq_along(plot.features)){
  plots[[i]] = VlnPlot(scRNA, group.by = "orig.ident", pt.size = 0,
                       features = plot.features[i]) + NoLegend()
}
violin <- wrap_plots(plots = plots, nrow = 2)
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 10, height = 8)


# rPCA integration 
scRNA.list <- SplitObject(scRNA, split.by="orig.ident")
scRNA.list <- lapply(X = scRNA.list, FUN = function(x) {
  # Normalization
  x <- NormalizeData(x, scale.factor = 10000)
  # Feature selection
  x <- FindVariableFeatures(x, x.low.cutoff=0.0125, x.high.cutoff=3,y.cutoff=0.5)
  # Scaling the data
  x <- ScaleData(x)
})


features <- SelectIntegrationFeatures(object.list = scRNA.list)
scRNA.list <- lapply(X = scRNA.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# RPCA integration regress out donor effects
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNA.list, anchor.features = features, reduction = "rpca")
scRNA.combined <- IntegrateData(anchorset = scRNA.anchors)

# Run the standard workflow for visualization and clustering
scRNA.combined <- ScaleData(scRNA.combined, verbose = FALSE)
scRNA.combined <- RunPCA(scRNA.combined, npcs = 30, verbose = FALSE)
scRNA.combined <- RunUMAP(scRNA.combined, dims = 1:30)

# Clustering
scRNA.combined <- FindNeighbors(scRNA.combined, dims = 1:5)
scRNA.combined <- FindClusters(scRNA.combined, resolution = 0.3)

# Annotation
new.cluster.ids <- c("Microglia","Astrocyte", "Neuron","Oligodendrocyte", "OPC", "Endothelial", "Unidentified", "Hybrid")
names(new.cluster.ids) <- levels(scRNA.combined)
scRNA.combined <- RenameIdents(scRNA.combined, new.cluster.ids)

# Visualization:
# APOE gene expression
dir.create("Visualization")
APOE <- FeaturePlot(scRNA.combined, features = "APOE")
ggsave("Visualization/APOE_Exp.pdf", plot = APOE, width = 10, height = 8)

# UMAP plot
UMAP <- DimPlot(scRNA.combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave("Visualization/UMAP_dimplot.pdf", plot = UMAP, width = 10, height = 8)

