process UTILS_CONVERT_BETWEEN_RDS_AND_H5AD {
    tag "convert_${file_basename}"
    label 'retry'
    label 'container_seuratdisk'

    input:
        path(file_to_convert)
    
    output:
        // Needs to return rds file for computing metrics
        path("${file_basename}.${converted_file_extension}")

    script:
        file_basename = file(file_to_convert).getSimpleName()
        file_extension = file(file_to_convert).extension
        
        converted_file_extension = file_extension == "rds" ? "h5ad" : "rds"

        """
        utils_convert_between_rds_and_h5ad.R \
        --input_path $file_to_convert
        """
}