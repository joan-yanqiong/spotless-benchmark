process COMPUTE_METRICS {
    tag "${method_name}_${output_suffix}"
    publishDir { "${params.outdir.metrics}/${output_suffix.replaceFirst(/_[a-z]{3}[0-9]+/, "")}" },
                mode: 'copy'
    debug true

    input:
        tuple val (method_name), path (props_file), path (sp_input)

    output:
        tuple val (method_name), path ("$metrics_file")

    script:
        output_suffix = file(sp_input).getSimpleName()
        metrics_file = "metrics_${method_name}_${output_suffix}${params.runID_metrics}"
        annot_remap = ( params.remap_annot ? "--remap $params.remap_annot" : "")

        """
        compute_metrics.R \
        --props_file $props_file --sp_input $sp_input \
        --output $metrics_file $annot_remap
        echo $metrics_file
        """
}
