process DSTG_CONVERT_REF_SEURAT_TO_MATRIX {
    label "retry"
    label "container_sp_dstg"

    input:
        path (sc_input)

    output:
        tuple path("${sc_filename}_matrix.rds"), path("${sc_filename}_label.rds"), emit: "rds"

    script:
        sc_filename = file(sc_input).getSimpleName()
        
        """
        dstg_convert_ref_seurat_to_matrix.R \
            --sc_input $sc_input --annot $params.annot
        """

}
