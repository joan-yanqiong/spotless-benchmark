#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(Seurat)
library(xbioc, quietly = TRUE)
library(Biobase)
library(magrittr)

# Replace default values by user input
par <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
# Convert numbers to numeric type
par[grepl("^[0-9\\.]+$", par)] <- as.numeric(par[grepl("^[0-9\\.]+$", par)])
print(par)

## START ##
cat("Reading input scRNA-seq reference from", par$sc_input, "\n")
seurat_obj_scRNA <- readRDS(par$sc_input)
DefaultAssay(seurat_obj_scRNA) <- "RNA"

if (prod(dim(seurat_obj_scRNA)) > 2**31) {
    cat("Reference is too large.\n")
    file.create("scrnaseq_ref_obj_downsampled.rds")
}
