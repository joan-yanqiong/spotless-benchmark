process MUSIC_CHECK_SCRNASEQ_REF_SIZE {
    label "mem_32G"
    label "time_1h"

    input:
        path sc_input

    output:
        path ("scrnaseq_ref_obj_downsampled.rds"), emit: "rds", optional: true
    
    script:
        args = (params.deconv_args.music ? params.deconv_args.music : "")
        """
        music_check_scrnaseq_ref_size.R \
            --sc_input $sc_input
        """
        
}
