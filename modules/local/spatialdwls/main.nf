process SPATIAL_DWLS_RUN {
    tag "spatialdwls_$sample_id"
    label "retry"
    label "time_2h"
    label "mem_32G"

    input:
        path(sc_input)
        tuple val(sample_id), path(sp_input)

    output:
        tuple val('spatialdwls'), val(sample_id), path("$output"), path(sp_input),  emit: "txt"

    script:
        output = "spatialdwls__${sample_id}.txt"

        args = (params.deconv_args.spatialdwls ? params.deconv_args.spatialdwls : "")
        """
        spatialdwls_run.R \
            --sc_input $sc_input \
            --sp_input $sp_input \
            --rootdir $projectDir \
            --annot $params.annot \
            --output $output $args
        """

}
