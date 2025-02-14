process CELL2LOCATION_FIT_MODEL {
    tag "c2l_${sample_id}"
    label "retry"
    label "container_cell2location"
    label "time_4h"
    label "mem_32G"
    label ( params.gpu ? "use_gpu" : "use_cpu" )
    debug true

    input:
        path (model)
        tuple val(sample_id), path(sp_input)

    output:
        tuple val('cell2location'), val(sample_id), path("$output")
    script:
        output = "cell2location__${sample_id}.preformat"
        epochs = ( params.epoch_fit ==~ /default/ ? "" : "-e $params.epoch_fit")
        args = ( params.deconv_args.cell2location.fit ? params.deconv_args.cell2location.fit : "" )
        cuda_device = ( params.gpu ? params.cuda_device : "cpu" )
        println ("Fitting cell2location model from file ${model} with ${ (params.gpu) ? "GPU" : "CPU" }...")
        println ("Arguments: ${args}")
        """
        source activate cell2loc_env
        cell2location_fit_model.py \
            $sp_input \
            $model \
            $cuda_device \
            $epochs \
            $args \
            -o \$PWD 
        mv proportions.tsv $output
        """
}
