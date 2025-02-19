library(Seurat)
library(SeuratObject)
library(magrittr)

##### SCRIPT TO CONVERT SEURAT OBJECT TO LIST WITH MATRIX FILE #####
par <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

cat("Reading input spatial data from", par$sp_input, "\n")
spatial_data <- readRDS(par$sp_input)

if (class(spatial_data) != "Seurat") {
    matrix_visium <- spatial_data$counts
} else { # If it is Seurat object
    DefaultAssay(spatial_data) <- names(spatial_data@assays)[grep("RNA|Spatial", names(spatial_data@assays))[1]]
    matrix_visium <- GetAssayData(spatial_data, slot = "counts")
}

colnames(matrix_visium) <- colnames(matrix_visium) %>%
    # Replace prohibited characters
    stringr::str_replace_all("[/ .\\-&]", "") %>%
    # Add "spot_" if colname starts with number (otherwise, error downstream)
    ifelse(stringr::str_detect(., "^[0-9]"),
        paste0("spot_", .), .
    )

cat("Saving new file as", paste0(par$sample_id, "_matrix.rds..."), "\n")
saveRDS(matrix_visium, paste0(par$sample_id, "_matrix.rds"))
