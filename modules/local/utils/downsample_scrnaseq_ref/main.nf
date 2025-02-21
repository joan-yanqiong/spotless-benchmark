process UTILS_DOWNSAMPLE_SCRNASEQ_REF {
    tag "convert_scrnaseq"
    label 'retry'
    label 'container_seuratdisk'

    input:
    path sc_input

    output:
    path "${output}", emit: rds, optional: true

    script:
    output = "scrnaseq_ref_obj_downsampled.rds"
    """
    utils_downsample_matrix.R \
        --sc_input ${sc_input} \
        --annot ${params.annot}
    """
}
