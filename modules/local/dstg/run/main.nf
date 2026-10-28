process DSTG_RUN {
    tag "dstg_${meta.sample_id}"
    label "retry"
    label "container_sp_dstg"

    input:
    tuple path(sc_input_matrix), path(sc_input_labels)
    tuple val(meta), path(sp_input)

    output:
    tuple val('dstg'), val(meta), path("${output}"), emit: "csv"

    script:
    output = "dstg__${meta.sample_id}.csv"

    """
        WORKDIR=\$PWD
        
        source activate dstg
        # Script doesn't run unless we're inside the directory
        cp -R ${params.deconv_args.dstg.dir}/DSTG/DSTG .
        cd DSTG/
        rm synthetic_data/*

        # Input: sc file, spatial file, and labels file
        Rscript convert_data.R ../${sc_input_matrix.getName} \
        ../${sp_input.getName} ../${sc_input_labels.getName}

        # Train the model
        python train.py

        # Add line numbers to output
        awk -v ln=1 '{print ln++  ","  \$0 }' DSTG_Result/predict_output.csv > tmp && mv tmp DSTG_Result/predict_output.csv
        
        # Add cell type labels
        { head -n 1 Datadir/Real_Label2.csv ; cat DSTG_Result/predict_output.csv; } > tmp && mv tmp DSTG_Result/predict_output.csv
        
        # Convert csv to tsv and move file to work directory
        sed -i -E 's/("([^"]*)")?,/\\2\t/g' DSTG_Result/predict_output.csv

        mv DSTG_Result/predict_output.csv ../${output}
    
        """

    stub:

    """
    touch ${output}
    """
}
