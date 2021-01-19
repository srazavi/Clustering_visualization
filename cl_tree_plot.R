library(clustree)
library(Seurat)
library(SeuratDisk)
library(tidyverse)

nn_suffix <- "_nn"
nn_number_vector <- c(10, 12)
clustering_res_prefix <- "lr_"

Convert("write/results.h5ad", dest = "h5seurat", overwrite = T)
scanpy_object_Seurat <- LoadH5Seurat(paste("write/results.h5seurat", sep=""))

dir.create("plots")
for (val in nn_number_vector) {
  nn_col_suffix <- paste("_", val, nn_suffix, sep="")
  clusters <- scanpy_object_Seurat@meta.data %>% select(ends_with(nn_col_suffix))
  colnames(clusters)
  clustree_graph <- 
    clustree(clusters, prefix=clustering_res_prefix, 
             suffix=nn_col_suffix, label_nodes=T)
  ggsave(paste("plots/clustering_tree_", nn_col_suffix, ".png", sep=""), clustree_graph)
}
