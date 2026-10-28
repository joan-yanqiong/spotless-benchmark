process STRIDE_BUILD_MODEL {
    tag "stride_build_$tag_suffix"
    label "retry"
    label "longer_time"
    label "container_stride"
    debug true

    input:
        path (sc_input)
    output:
        tuple path ("model.zip"), path("txt_files.zip")

    script:
        args = ( params.deconv_args.stride ? params.deconv_args.stride : "" )
        tag_suffix = file(sc_input).getSimpleName()
        
        """
        source activate stride
        echo "Creating annotation and dummy spatial files..."
        python $projectDir/subworkflows/deconvolution/stride/createFiles.py \
           $sc_input $params.annot dummy_arg

        STRIDE deconvolve --sc-count $sc_input --sc-celltype annot.txt \
            --st-count dummy_st.h5ad $args

        # The model/ folder contains multiple files, so zip it
        zip -r model.zip model
        zip txt_files.zip Model_selection.txt Gene_dict.txt
        """
}

