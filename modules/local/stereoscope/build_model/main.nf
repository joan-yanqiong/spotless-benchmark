process STEREOSCOPE_BUILD_MODEL {
    tag "stereo_build_$tag_suffix"
    label "retry"
    label "longer_time"
    label ( params.gpu ? "use_gpu" : "use_cpu" )
    label "container_stereoscope"
    debug true

    input:
        path (sc_input)
    output:
        tuple path ("R*.tsv"), path ("logits*.tsv")
    // TODO handle compute_hvgs and build_model
    script:
        tag_suffix = file(sc_input).getSimpleName()
        epochs = ( params.epoch_build ==~ /default/ ? "" : "-sce $params.epoch_build")
        args = ( params.deconv_args.stereoscope ? params.deconv_args.stereoscope : "" )
        gpu_flag = ( params.gpu ? "--gpu" : "" )
        println ("Building stereoscope model with ${ (params.gpu) ? "GPU" : "CPU" }...")

        """
        source activate stereoscope
        export CUDA_VISIBLE_DEVICES=$params.cuda_device

        echo "Arguments: $args"
        if [[ "$args" == *"-gl"* ]] && [[ "$args" == *"-hvg"* ]]
        then
          echo "ERROR: Both gene list and hvg option is provided. You can only choose one."
          exit 1

        elif [[ "$args" == *"-hvg"* ]]
        then
          echo "Computing HVGs..."
          stereoscope_compute_hvgs.py \
            --sc_cnt $sc_input $args
          
          # Remove -hvg from argument, add -gl
          re='(.*)-hvg[ ]*[0-9]+(.*)'
          [[ "$args" =~ \$re ]]
          new_args="\${BASH_REMATCH[1]}\${BASH_REMATCH[2]} -gl hvgs.txt"

        else
          new_args="$args"
        fi

        stereoscope run --sc_cnt $sc_input --label_colname $params.annot \
            $epochs \$new_args $gpu_flag -o \$PWD
        """
}
