process SEURAT_RUN {
    tag "seurat_$output_suffix"

    publishDir { "${params.outdir.props}/${output_suffix.replaceFirst(/_[a-z]{3}[0-9]+/, "")}" },
                mode: 'copy', pattern: "proportions_*"

    input:
        tuple path (sc_input), path (sp_input)
    output:
        tuple val('seurat'), path("$output"), path (sp_input)
    
    script:
        output_suffix = file(sp_input).getSimpleName()
        output = "proportions_seurat_${output_suffix}${params.runID_props}"
        args = (params.deconv_args.seurat ? params.deconv_args.seurat : "")

        """
        seurat_run.R \
            --sc_input $sc_input --sp_input $sp_input \
            --annot $params.annot --output $output $args
        """

}
