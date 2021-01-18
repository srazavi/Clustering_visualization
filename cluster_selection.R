library(clustree)
library(Seurat)
library(SeuratDisk)
library(tidyverse)

nn_suffix <- "_nn"
nn_number_vector <- c(10, 12)
clustering_res_prefix <- "lr_"

Convert("write/results.h5ad", dest = "h5seurat", overwrite = T)
scanpy_object_Seurat <- LoadH5Seurat(paste("write/results.h5seurat", sep=""))

nn_col_suffix <- "_10_nn"
clusters_10_nn <- scanpy_object_Seurat@meta.data %>% select(ends_with(nn_col_suffix))
colnames(clusters_10_nn)
clustree_graph <- 
  clustree(clusters_10_nn, prefix=clustering_res_prefix, suffix=nn_col_suffix, label_nodes=T)
clustree_graph
