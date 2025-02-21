#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(SeuratDisk)
library(Seurat)
options(Seurat.object.assay.version = "v3")

print(getOption("Seurat.object.assay.version"))

par <- list(
    output_file_ext = ""
)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
par[names(args)] <- args

file_name <- stringr::str_split(basename(par$input_path), "\\.")[[1]][1]
ext <- stringr::str_split(basename(par$input_path), "\\.")[[1]][2]

# Convert RDS to H5AD
if (tolower(ext) == "rds") {
    message("Input file is rds, will convert to h5ad...\n")
    message("Reading input ", par$input_path, "\n")
    input_obj <- readRDS(par$input_path)

    # Synthspot object
    if (class(input_obj) != "Seurat") {
        seurat_obj <- CreateSeuratObject(counts = input_obj$counts)
    } else {
        seurat_obj <- input_obj
    }
    rm(input_obj)

    # Use raw counts
    def_assay <- names(seurat_obj@assays)[grep(
        "RNA|Spatial",
        names(seurat_obj@assays)
    )[1]]
    message(paste0("Set default assay: ", def_assay, "..."))
    DefaultAssay(seurat_obj) <- def_assay

    seurat_obj_version <- as.character(seurat_obj@version)
    message(paste0(
        "Detected Seurat object version: ",
        seurat_obj_version,
        "...",
        sep = " "
    ))
    # SeuratDisk cannot work with a Seurat object older than v3.1.2
    if (compareVersion(seurat_obj_version, "3.1.2") == -1) {
        message("Seurat object is too old, creating a new one...")
        seurat_obj <- CreateSeuratObject(
            counts = GetAssayData(seurat_obj, slot = "counts"),
            assay = def_assay,
            meta.data = seurat_obj@meta.data
        )
        message("Seurat object created...")

        # If your Seurat object is a v5, then create a new object
    } else if (startsWith(seurat_obj_version, "5")) {
        print(seurat_obj)
        message("Create a new Seurat object...")
        seurat_obj <- CreateSeuratObject(
            counts = CreateAssayObject(counts = seurat_obj[[def_assay]]@counts),
            assay = def_assay,
            meta.data = seurat_obj@meta.data
        )
        message("Seurat object created...")
    }

    # If the object has been preprocessed before, SeuratDisk is going to place the
    # raw counts in the .raw.X instead of .X, which wouldn't work downstream
    # because the Python methods access .X instead of .raw.X
    # So we check if the scale matrix exists, and if the counts and data slots are different
    message("Check if the Seurat object has been preprocessed before...")
    if (
        !isTRUE(all.equal(
            GetAssayData(seurat_obj, slot = "counts"),
            GetAssayData(seurat_obj, slot = "data")
        )) ||
            all(dim(GetAssayData(seurat_obj, slot = "scale.data")) > 0)
    ) {
        message(
            "The Seurat object has been preprocessed. Creating a new object..."
        ) # nolint
        seurat_obj <- CreateSeuratObject(
            counts = GetAssayData(seurat_obj, slot = "counts"),
            assay = DefaultAssay(seurat_obj),
            meta.data = seurat_obj@meta.data
        )
    }
    gc()

    if (file.exists(paste0(file_name, ".h5ad"))) {
        return("h5ad file already exists")
    }

    message("Writing h5Seurat file...\n")
    SaveH5Seurat(seurat_obj, filename = paste0(file_name, ".h5seurat"))

    message("Writing", paste0(file_name, ".h5ad"), "\n")
    Convert(paste0(file_name, ".h5seurat"), dest = "h5ad")
    file.remove(paste0(file_name, ".h5seurat"))

    # Convert H5AD to RDS
} else if (tolower(ext) == "h5ad") {
    message("Input file is h5ad, will convert to rds...\n")
    dirname(Sys.readlink(par$input_path))
    print(par$input_path)
    message("Writing h5Seurat file...\n")

    # h5seurat is saved at original directory (follow symlink)
    h5seurat_path <- file.path(
        dirname(Sys.readlink(par$input_path)),
        paste0(file_name, ".h5seurat")
    )

    if (file.exists(h5seurat_path)) {
        message("h5seurat file already exists, using the existing one...\n")
    } else {
        Convert(par$input_path, dest = "h5seurat")
    }

    message("Writing", paste0(file_name, par$output_file_ext, ".rds"), "\n")
    h5seurat_file <- LoadH5Seurat(h5seurat_path)
    saveRDS(h5seurat_file, paste0(file_name, par$output_file_ext, ".rds"))
    file.remove(h5seurat_path)
}

message("Finished!")
