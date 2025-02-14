process UTILS_FORMAT_TSV_FILE {
    tag "format_${method_name}"
    label 'trivial'
    label 'container_tidyverse'
    label 'time_30m'
    label 'mem_16G'
    input:
        tuple val(method_name), val(sample_id), path (tsv_file)
    output:
        tuple val(method_name), val(sample_id), path (new_tsv_file)
    script:
        new_tsv_file = "${method_name}__${sample_id}.txt"
        template 'utils_format_tsv_file.R'

}
