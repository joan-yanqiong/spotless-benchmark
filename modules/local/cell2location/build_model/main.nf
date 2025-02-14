process CELL2LOCATION_BUILD_MODEL {
    tag "c2l_build_$tag_suffix"
    label "retry"
    label "time_4h"
    label "mem_32G"
    label "container_cell2location"
    label ( params.gpu ? "use_gpu" : "use_cpu" )
    debug true

    input:
        path (sc_input)
    output:
        path "sc.h5ad", emit: h5ad

    script:
        tag_suffix = file(sc_input).getSimpleName()
        sample_id_arg = ( params.sampleID ==~ /none/ ? "" : "-s $params.sampleID" )
        epochs = ( params.epoch_build ==~ /default/ ? "" : "-e $params.epoch_build")
        args = ( params.deconv_args.cell2location.build ? params.deconv_args.cell2location.build : "" )
        cuda_device = ( params.gpu ? params.cuda_device : "cpu" )
        println ("Building cell2location model with ${ (params.gpu) ? "GPU" : "CPU" }...")
        """
        source activate cell2loc_env
        cell2location_build_model.py \
            $sc_input $cuda_device -a $params.annot $sample_id_arg $epochs $args -o \$PWD 
        """

}
