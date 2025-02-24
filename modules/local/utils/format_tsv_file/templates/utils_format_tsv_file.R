#!/usr/bin/env Rscript
deconv_matrix <- read.table(
    "$tsv_file",
    sep = "\t",
    header = TRUE,
    row.names = 1
)
colnames(deconv_matrix) <- stringr::str_replace_all(
    colnames(deconv_matrix),
    "[/ .]",
    ""
)
deconv_matrix <- deconv_matrix[, sort(
    colnames(deconv_matrix),
    method = "shell"
)]
write.table(
    deconv_matrix,
    file = "$new_tsv_file",
    sep = "\t",
    quote = FALSE,
    row.names = TRUE
)
