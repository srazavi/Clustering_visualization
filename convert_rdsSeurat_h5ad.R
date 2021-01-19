library(Seurat)
library(SeuratDisk)

data <- readRDS(file = "org.combinedR2_clusters_names.rds")
# data@assays$RNA@counts@Dimnames$

SaveH5Seurat(data, filename = "org.combinedR2_clusters_names.h5Seurat")
Convert("org.combinedR2_clusters_names.h5Seurat", dest = "h5ad", overwrite = T)

