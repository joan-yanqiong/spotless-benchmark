process UTILS_CREATE_DUMMY_FILE {
    tag "dummy_${file_name}"
    label 'trivial'

    input:
        path (h5ad_file)
    output:
        tuple path (h5ad_file), path ("${file_name}.rds")

    script:
        file_name = file(h5ad_file).getSimpleName()
        """
        touch ${file_name}.rds
        """
}