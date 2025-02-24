#!/usr/bin/env Rscript

## HELPER FUNCTIONS TO DOWNSAMPLE REFERENCE SCRNA-SEQ OBJECT ##
# Used in MuSiC and SpatialDWLS
# When the reference exceeds 2**31 elements, cells and genes
# will be downsampled to not exceed the dense matrix limit

Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(Seurat)
library(magrittr)

par <- list(
    # Downsampling cells
    downsample_cells = TRUE, # If dense matrix is too big, downsample cells
    target_n_cells = 10000, # Max no. of cells per cell type

    # Downsampling genes
    downsample_genes = TRUE, # If dense matrix is too big, downsample genes
    n_hvgs = 3000, # Number of highly variable genes to keep
    pct = 0.1, # Percentage of cells which genes have to be expressed
    assay_oi = "RNA",

    # Filtering spatial object
    filter_spots = "none"
)

# Replace default values by user input
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
par[names(args)] <- args
# Convert numbers to numeric type
par[grepl("^[0-9\\.]+$", par)] <- as.numeric(par[grepl("^[0-9\\.]+$", par)])
print(par)

# This function returns a vector of genes that are expressed in at least
# pct fraction of the celltype of interest

#' @title Get expressed genes
#' @param celltype_oi celltype of interest
#' @param seurat_obj seurat object containing the cell type of interest
#' @param assay_oi assay to use (default='RNA')
#' @param pct pct (fraction; value between 0 and 1) of the cells that need to express a gene to be kept (default=0.1)
#' @return vector of genes that are expressed in at least pct fraction of the cell type of interest
get_expressed_genes <- function(
    celltype_oi,
    seurat_obj,
    annot_col,
    assay_oi = "RNA",
    pct = 0.1
) {
    # Subset expression matrix to contain only celltype of interest
    cells_oi <- colnames(seurat_obj)[seurat_obj[[annot_col]] == celltype_oi]
    exprs_mat <- seurat_obj[[assay_oi]]@data %>% .[, cells_oi]

    n_cells_oi <- ncol(exprs_mat)
    if (n_cells_oi < 5000) {
        genes <- exprs_mat %>%
            apply(1, function(x) {
                sum(x > 0) / n_cells_oi
            }) %>%
            .[. >= pct] %>%
            names()
    } else {
        # If there are more than 5000 cells there seems to be some memory issue
        # Split into chunks of 100 genes
        nchunk_size <- 100
        splits <- split(
            seq_len(nrow(exprs_mat)),
            ceiling(seq_along(seq_len(nrow(exprs_mat))) / nchunk_size)
        )
        genes <- splits %>%
            lapply(
                function(genes_indices, exprs, pct, n_cells_oi) {
                    begin_i <- genes_indices[1]
                    end_i <- genes_indices[length(genes_indices)]
                    exprs <- exprs[begin_i:end_i, ]
                    genes <- exprs %>%
                        apply(1, function(x) {
                            sum(x > 0) / n_cells_oi
                        }) %>%
                        .[. >= pct] %>%
                        names()
                },
                exprs_mat,
                pct,
                n_cells_oi
            ) %>%
            unlist() %>%
            unname()
    }
    return(genes)
}


#' @title Downsampling cells
#' @description This function returns a list of cell indices such that each cell type contains a maximum of target_n_cells
#' @param seurat_obj Seurat object
#' @param annot_col column in metadata of Seurat objct with annotation
#' @param target_n_cells max. no. of cells targeted to be sampled for each cell type (default=10000)
#' @return list of indices
get_downsampled_cells <- function(
    seurat_obj,
    annot_col,
    target_n_cells = 10000
) {
    index_keep <- sapply(
        unique(seurat_obj[[annot_col, drop = TRUE]]),
        function(celltype) {
            indices_oi <- which(seurat_obj[[annot_col]] == celltype)
            n_cells <- min(target_n_cells, length(indices_oi))
            sample(indices_oi, n_cells, replace = FALSE)
        }
    )
    return(unlist(index_keep) %>% sort())
}

## START ##
message("Reading input scRNA-seq reference from", par$sc_input, "\n")
seurat_obj_scRNA <- readRDS(par$sc_input)
DefaultAssay(seurat_obj_scRNA) <- "RNA"
ncelltypes <- length(unique(seurat_obj_scRNA[[par$annot, drop = TRUE]]))
message("Found ", ncelltypes, "cell types in the reference.\n")

if (prod(dim(seurat_obj_scRNA)) <= 2**31) {
    message("Reference size is OK...")
    file.create("ok.rds")
} else {
    message("Reference is too large, downsample...")
    ### MAIN CODE ###
    new_cells <- colnames(seurat_obj_scRNA)
    features_keep <- rownames(seurat_obj_scRNA)

    if (par$downsample_cells) {
        message(paste(
            "Downsampling to",
            par$target_n_cells,
            "cells per cell type.."
        ))
        new_cells <- get_downsampled_cells(
            seurat_obj_scRNA,
            par$annot,
            par$target_n_cells
        )
        message(paste("Reference now has", length(new_cells), "cells..."))
    }

    if (par$downsample_genes) {
        message(paste("Downsampling genes by keeping", par$n_hvgs, "HVGS..."))
        message(paste(
            "and genes that are expressed in ",
            par$pct,
            " fraction per cell type..."
        ))

        # Preprocess reference object to get HVGs and log normalized data
        seurat_obj_scRNA <- seurat_obj_scRNA %>%
            NormalizeData() %>%
            FindVariableFeatures(nfeatures = par$n_hvgs)

        var_genes <- VariableFeatures(seurat_obj_scRNA)
        expressed_genes_list <- unique(seurat_obj_scRNA[[
            par$annot,
            drop = TRUE
        ]]) %>%
            lapply(., function(celltype) {
                get_expressed_genes(
                    celltype,
                    seurat_obj_scRNA,
                    par$annot,
                    par$assay_oi,
                    par$pct
                )
            })

        features_keep <- union(
            var_genes,
            expressed_genes_list %>%
                unlist() %>%
                unique()
        )
        message("Reference now has", length(features_keep), "genes.\n")
    }

    # Use as.numeric to prevent integer overflow (`length` returns integer)
    if (
        as.numeric(length(features_keep)) * as.numeric(length(new_cells)) >
            2**31
    ) {
        stop(paste0(
            "Reference is still larger than 2^31 with ",
            as.numeric(length(features_keep)),
            "x",
            as.numeric(length(new_cells)),
            " elements..."
        ))
    }

    seurat_obj_scRNA <- seurat_obj_scRNA[features_keep, new_cells]
    message("Save downsampled Seurat object as RDS...")
    saveRDS(seurat_obj_scRNA, "scrnaseq_ref_obj_downsampled.rds")
    gc()
}
message("Finished!")
