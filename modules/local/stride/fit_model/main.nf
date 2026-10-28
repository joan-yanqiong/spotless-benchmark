
process STRIDE_FIT_MODEL {
    tag "stride_$sp_file_basename"
    label "retry"
    label "container_stride"
    debug true

    input:
        tuple path (sp_input), path (sp_input_rds)
        tuple path (model_zip), path(txt_files_zip)
    output:
        tuple val('stride'), path("$output"), path (sp_input_rds)
    script:
        sp_file_basename = file(sp_input).getSimpleName()
        output = "proportions_stride_${sp_file_basename}${params.runID_props}.preformat"
        args = ( params.deconv_args.stride ? params.deconv_args.stride : "" )
        
        """
        source activate stride
        unzip $model_zip
        unzip $txt_files_zip
        
        STRIDE deconvolve --model-dir model/ --st-count $sp_input $args
        
        mv STRIDE_spot_celltype_frac.txt $output
        """
}

