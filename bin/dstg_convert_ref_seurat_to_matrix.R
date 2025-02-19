library(Seurat)
library(SeuratObject)
library(magrittr)

##### SCRIPT TO CONVERT SEURAT OBJECT TO LIST WITH MATRIX FILE #####
par <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

cat("Reading input scRNA-seq reference from", par$sc_input, "\n")
seurat_obj_scRNA <- readRDS(par$sc_input)
DefaultAssay(seurat_obj_scRNA) <- "RNA"
matrix_scRNA <- GetAssayData(seurat_obj_scRNA, slot = "counts")

file_name_scrna <- stringr::str_split(basename(par$sc_input), "\\.")[[1]][1]
cat("Saving new file as", paste0(file_name_scrna, "_matrix.rds...\n"))
saveRDS(matrix_scRNA, paste0(file_name_scrna, "_matrix.rds"))

ncelltypes <- length(unique(seurat_obj_scRNA[[par$annot, drop = TRUE]]))
cat("Found ", ncelltypes, "cell types in the reference.\n")
cell_types <- stringr::str_replace_all(
    seurat_obj_scRNA[[par$annot, drop = TRUE]],
    "[/ .\\-&]", ""
) # Replace prohibited characters
names(cell_types) <- colnames(seurat_obj_scRNA)
cat("Saving cell types as", paste0(file_name_scrna, "_label.rds..."), "\n")
saveRDS(cell_types, paste0(file_name_scrna, "_label.rds"))
