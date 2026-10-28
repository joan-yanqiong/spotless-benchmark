process NNLS_RUN {
    tag "nnls_$output_suffix"
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2
    publishDir { "${params.outdir.props}/${output_suffix.replaceFirst(/_[a-z]{3}[0-9]+/, "")}" },
                mode: 'copy', pattern: "proportions_*"

    input:
        tuple path (sc_input), path (sp_input)
    output:
        tuple val('nnls'), path("$output"), path (sp_input)
    
    script:
        output_suffix = file(sp_input).getSimpleName()
        output = "proportions_nnls_${output_suffix}${params.runID_props}"

        """
        nnls_run.R \
            --sc_input $sc_input --sp_input $sp_input \xw
            --annot $params.annot --output $output
        """

}
