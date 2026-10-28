process SPOTLIGHT_RUN {
    tag "spotlight_$sample_id"
    label "longer_time"
    label "retry"
    label "mem_32G"
    label "time_1h"

    input:
        path (sc_input)
        tuple val(sample_id), path (sp_input)

    output:
        tuple val('spotlight'), val(sample_id), path("$output")

    script:
        output = "spotlight__${sample_id}"
        args = (params.deconv_args.spotlight ? params.deconv_args.spotlight : "")
        """
        spotlight_run.R \
            --sc_input $sc_input --sp_input $sp_input \
            --annot $params.annot --output $output $args
        """

}