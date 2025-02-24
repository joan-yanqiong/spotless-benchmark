process MUSIC_DOWNSAMPLE_SCRNASEQ_REF {
    label "mem_32G"
    label "time_1h"

    input:
        path sc_input

    when: sc_input.size() > 0

    output:
        path ("scrnaseq_ref_obj_downsampled.rds"), emit: "rds", optional: true
    
    script:
        args = (params.deconv_args.music ? params.deconv_args.music : "")
        """
        utils_downsample_matrix.R \
            --sc_input $sc_input \
            --annot $params.annot $args

        """
}
