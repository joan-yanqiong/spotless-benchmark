process STEREOSCOPE_FIT_MODEL {
    tag "stereo_$sp_file_basename"
    label "retry"
    label "longer_time"
    label ( params.gpu ? "use_gpu" : "use_cpu" )
    label "container_stereoscope"
    debug true

    input:
        tuple path (sp_input), path (sp_input_rds)
        tuple path (r_file), path (logits_file)
    output:
        tuple val('stereoscope'), path("$output"), path (sp_input_rds)
    script:
        sp_file_basename = file(sp_input).getSimpleName()
        output = "proportions_stereoscope_${sp_file_basename}${params.runID_props}.preformat"
        epochs = ( params.epoch_fit ==~ /default/ ? "" : "-ste $params.epoch_fit")
        args = ( params.deconv_args.stereoscope ? params.deconv_args.stereoscope : "" )
        gpu_flag = ( params.gpu ? "--gpu" : "" )
        args = args.replaceFirst(/-hvg[ ]*[0-9]+/, "") // Replace HVG in original argument

        println ("Received model files $r_file and $logits_file")
        println ("Fitting stereoscope model with ${ (params.gpu) ? "GPU" : "CPU" }...")
        
        """
        source activate stereoscope
        export CUDA_VISIBLE_DEVICES=$params.cuda_device
        
        stereoscope run --sc_fit $r_file $logits_file \
            --st_cnt $sp_input $epochs $args $gpu_flag -o \$PWD
        mv $sp_file_basename/W*.tsv $output
        """
}
