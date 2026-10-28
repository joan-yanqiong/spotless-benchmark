process MUSIC_RUN {
    tag "music_$sample_id"
    label "mem_32G"
    label "time_4h"

    input:
        path sc_input
        tuple val(sample_id), path (sp_input)
    output:
        tuple val('music'), val(sample_id), path("$output"), emit: "txt"
    
    script:
        output = "music__${sample_id}"
        args = (params.deconv_args.music ? params.deconv_args.music : "")

        """
        music_run.R \
            --sc_input $sc_input --sp_input $sp_input \
            --annot $params.annot --output $output \
            --rootdir $projectDir \
            --sampleID $params.sampleID $args
        """

}
