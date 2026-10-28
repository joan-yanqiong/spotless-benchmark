process DESTVI_BUILD_MODEL {
    tag "destvi_build_$tag_suffix"
    label "retry"
    label "longer_time"
    label ( params.gpu ? "use_gpu" : "use_cpu" )
    label "container_destvi"
    debug true

    input:
        path (sc_input)
        tuple path (sp_input_h5ad), path (sp_input_rds)
    output:
        tuple path ("adata.h5ad"), path ("model.pt")

    script:
        tag_suffix = file(sc_input).getSimpleName()
        epochs = ( params.epoch_build ==~ /default/ ? "" : "-e $params.epoch_build")
        args = ( params.deconv_args.destvi.build ? params.deconv_args.destvi.build : "" )
        cuda_device = ( params.gpu ? params.cuda_device : "cpu" )
        println ("Building DestVI model with ${ (params.gpu) ? "GPU" : "CPU" }...")
        """
        source activate scvi-env
        build_model.py \
            $sc_input $sp_input_h5ad $cuda_device -a $params.annot \
            -o \$PWD/model $epochs $args
        mv model/* .
        """

}
