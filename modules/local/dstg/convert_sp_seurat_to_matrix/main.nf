process DSTG_CONVERT_SP_SEURAT_TO_MATRIX {
    tag "dstg_${meta.sample_id}"
    label "retry"
    label "container_sp_dstg"

    input:
    tuple val(meta), path(sp_input)

    output:
    tuple val(meta), path(output), emit: "rds"

    script:
    output = "${meta.sample_id}_matrix.rds"
    """
        dstg_convert_sp_seurat_to_matrix.R \
            --sp_input ${sp_input} \
            --sample_id ${meta.sample_id}
        """

    stub:
    """
    touch ${output}
    """
}
