process RCTD_RUN {
    tag "rctd_$sample_id"
    label "retry"
    label "time_2h"
    label "mem_32G"

    input:
        path(sc_input)
        tuple val(sample_id), path(sp_input)

    output:
        tuple  val('rctd'), val(sample_id), path("$output"), emit: "txt"

    script:
        output = "rctd__${sample_id}.txt"
        num_cores = task.cpus
        args = (params.deconv_args.rctd ? params.deconv_args.rctd : "")

        """
        rctd_run.R \
            --sc_input $sc_input \
            --sp_input $sp_input \
            --annot $params.annot \
            --output $output \
            --num_cores ${num_cores} $args
        """

}
