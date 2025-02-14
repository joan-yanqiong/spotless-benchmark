// This is for if you want to convert files outside the normal workflow
process UTILS_CONVERT_DATA_FORMAT_FOR_WORKFLOW {
    tag "convert_${file_basename}"
    container 'csangara/seuratdisk:latest'
    publishDir { "${file_directory}/" },
               mode: 'copy', pattern: "*.${ext_to_include}"
    label 'retry'
    label 'container_seuratdisk'

    input:
        path (file_to_convert)
    
    output:
        // Needs to return rds file for computing metrics
        tuple path ("${file_basename}.h5ad"), path ("${file_basename}${output_file_ext}.rds")

    script:
        file_basename = file(file_to_convert).getSimpleName()
        file_directory = params.convert_input =~ /\*/ ? file(params.convert_input)[0].getParent() :
                                                        file(params.convert_input).getParent()
        ext_to_include = file(file_to_convert).getExtension().toLowerCase() =~ /h5/ ? "rds" : "h5ad"
        output_file_ext = ext_to_include ==~ /rds/ ? "_fromH5AD" : ""
        println("Saving ${ext_to_include} file to ${file_directory}")
        """
        utils_convert_between_rds_and_h5ad.R \
        --input_path $file_to_convert --output_file_ext $output_file_ext
        """
}
