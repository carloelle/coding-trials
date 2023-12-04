library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)



image<-Read10X_Image('/Users/leonardi.carlo/Desktop/A1_outs/Spatial1/outs/spatial')
SpatialA1<-Load10X_Spatial('/Users/leonardi.carlo/Desktop/A1_outs/Spatial1/outs')

plot1 <- VlnPlot(SpatialA1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(SpatialA1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

SpatialA1 <- SCTransform(SpatialA1, assay = "Spatial", verbose = FALSE)

{# rerun normalization to store sctransform residuals for all genes
SpatialA1_1 <- SCTransform(SpatialA1_1, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
# also run standard log normalization for comparison
SpatialA1_1 <- NormalizeData(SpatialA1_1, verbose = FALSE, assay = "Spatial")

# Computes the correlation of the log normalized data and sctransform residuals with the
# number of UMIs
SpatialA1_1 <- GroupCorrelation(SpatialA1_1, group.assay = "Spatial", assay = "Spatial", slot = "data", do.plot = FALSE)
SpatialA1_1 <- GroupCorrelation(SpatialA1_1, group.assay = "Spatial", assay = "SCT", slot = "scale.data", do.plot = FALSE)

p1 <- GroupCorrelationPlot(SpatialA1_1, assay = "Spatial", cor = "nCount_Spatial_cor") + ggtitle("Log Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
p2 <- GroupCorrelationPlot(SpatialA1_1, assay = "SCT", cor = "nCount_Spatial_cor") + ggtitle("SCTransform Normalization") +
  theme(plot.title = element_text(hjust = 0.5))
p1 + p2}

SpatialA1 <- SCTransform(SpatialA1, assay = "Spatial", verbose = FALSE)
SpatialA1 <- RunPCA(SpatialA1, assay = "SCT", verbose = FALSE)
SpatialA1 <- FindNeighbors(SpatialA1, reduction = "pca", dims = 1:30)
SpatialA1 <- FindClusters(SpatialA1, verbose = FALSE)
SpatialA1 <- RunUMAP(SpatialA1, reduction = "pca", dims = 1:30)

p3 <- DimPlot(SpatialA1, reduction = "umap", label = TRUE)
p4 <- SpatialDimPlot(SpatialA1, label = TRUE, label.size = 3)
p3 + p4


#macro markers
SpatialFeaturePlot(SpatialA1, features = c("Adgre1",
'Cd163',
'Mrc1',
'Itgax',
'Cd63'))
#t cell
SpatialFeaturePlot(SpatialA1, features = c('Cd3e',
'Cd4',
'Cd8a',
'Ifng'))



SpatialA1 <- FindSpatiallyVariableFeatures(SpatialA1, assay = "SCT", features = VariableFeatures(SpatialA1)[1:1000],
                                       selection.method = "markvariogram")

top.features <- head(SpatiallyVariableFeatures(SpatialA1, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(SpatialA1, features = top.features, ncol = 3, alpha = c(0.1, 1))


