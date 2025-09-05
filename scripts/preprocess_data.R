# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("SeuratObject")
# install.packages("Seurat")
# install.packages("tidyverse")
# install.packages("patchwork")
# install.packages("bpCell")
# install.packages("hdf5r")

# install.packages("remotes")
# remotes::install_github("mojaveazure/seurat-disk")

library(SeuratDisk)
library(Seurat)
library(Matrix)
library(tools)

output_dir <- "converted_files_to_python_readable"
dir.create(output_dir, showWarnings = FALSE)

rds_files <- c(
  "data/TMA1_minimal.rds",
  "data/TMA2_minimal.rds",
  "data/TMA3_minimal.rds",
  "data/TMA4_minimal.rds"
)


process_rds <- function(file_path) {
  file_name <- file_path_sans_ext(basename(file_path))
  message("Processing: ", file_name)
  
  file_dir <- file.path(output_dir, file_name)
  dir.create(file_dir, showWarnings = FALSE)
  
  data <- readRDS(file_path)
  
  ## 1️⃣ Counts matrix
  message("Counts matrix: ", file_name)
  counts_mat <- GetAssayData(data, assay = "Nanostring", layer = "counts")
  writeMM(counts_mat, file = file.path(file_dir, "counts.mtx"))
  
  ## 2️⃣ Normalized data (optional, useful in Python)
  message("Normalized Data : ", file_name)
  if ("data" %in% slotNames(data@assays$Nanostring)) {
    norm_data <- GetAssayData(data, assay = "Nanostring", layer = "data")
    writeMM(norm_data, file = file.path(file_dir, "normalized_data.mtx"))
  }
  
  ## 3️⃣ Genes & cells
  message("Genes data: ", file_name)
  write.csv(rownames(counts_mat), file = file.path(file_dir, "genes.csv"), row.names = FALSE)
  message("Cell data: ", file_name)
  write.csv(colnames(counts_mat), file = file.path(file_dir, "cells.csv"), row.names = FALSE)
  
  ## 4️⃣ Metadata
  message("Meta data: ", file_name)
  write.csv(data@meta.data, file = file.path(file_dir, "meta_data.csv"))
  
  ## 5️⃣ Spatial coordinates
  message("Spatial Coordination data: ", file_name)
  if ("images" %in% slotNames(data) && length(data@images) > 0) {
    coords <- GetTissueCoordinates(data)
    write.csv(coords, file = file.path(file_dir, "spatial_coords.csv"))
  }
  
  ## 6️⃣ Cluster identities
  message("cluster data: ", file_name)
  write.csv(Idents(data), file = file.path(file_dir, "cluster_identities.csv"))
  
  ## 7️⃣ Dimensional reductions
  message("Dimensional reduction data: ", file_name)
  if (length(Reductions(data)) > 0) {
    for (red in names(Reductions(data))) {
      emb <- Embeddings(data, reduction = red)
      write.csv(emb, file = file.path(file_dir, paste0("reduction_", red, ".csv")))
    }
  }
  
  ## 8️⃣ Segmentation boundaries (if available)
  message("Segmentation boundary: ", file_name)
  if ("segmentation" %in% slotNames(data@images[[1]])) {
    seg <- data@images[[1]]@segmentation
    if (!is.null(seg)) {
      write.csv(seg@coordinates, file = file.path(file_dir, "segmentation_coords.csv"))
    }
  }
  
  ## 9️⃣ Molecule-level data (if available)
  message("Molecule-level data: ", file_name)
  mol_list <- data@images[[1]]@molecules[[1]]
  mol_df <- do.call(rbind, lapply(names(mol_list), function(cell) {
    df <- as.data.frame(mol_list[[cell]])
    df$cell <- cell
    df
  }))
  write.csv(mol_df, file.path(file_dir, "molecules.csv"), row.names = FALSE)
  
  message("✅ Saved in: ", file_dir)
}

# Apply function on every rds file
lapply(rds_files, process_rds)




#*****************************************************************************
# Read data for each TMA

# Load the .rds file
data <- readRDS("data/TMA1_minimal.rds")

counts_mat <- GetAssayData(data, assay = "Nanostring", layer = "counts")
write.csv(counts_mat, file = "gene_expression_tma1.csv")


writeMM(counts_mat, file = "counts.mtx")
write.csv(rownames(counts_mat), file = "genes.csv")
write.csv(colnames(counts_mat), file = "cells.csv")
write.csv(data@meta.data, "meta_data.csv")
coords <- GetTissueCoordinates(data, filtered=True)
write.csv(coords, "spatial_coords.csv")

centroids_cell_df <- as.data.frame(data@images[[1]]@boundaries[["centroids"]]@cells)
centroids_coords_df <- as.data.frame(data@images[[1]]@boundaries[["centroids"]]@coords)
merged_df <- cbind(centroids_cell_df, centroids_coords_df)
write.csv(merged_df, "centroids_merged.csv", row.names = FALSE)

#*****************************************************************************