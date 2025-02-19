process RCTD_RUN {
    tag "rctd_${meta.sample_id}"
    label "retry"
    label "time_2h"
    label "mem_32G"

    input:
    path sc_input
    tuple val(meta), path(sp_input)

    output:
    tuple val('rctd'), val(meta), path("${output}"), emit: "txt"

    script:
    output = "rctd__${meta.sample_id}.txt"
    num_cores = task.cpus
    args = (params.deconv_args.rctd ? params.deconv_args.rctd : "")

    """
        rctd_run.R \
            --sc_input ${sc_input} \
            --sp_input ${sp_input} \
            --annot ${params.annot} \
            --output ${output} \
            --num_cores ${num_cores} ${args}
        """

    stub:
    """
    touch ${output}
    """
}
