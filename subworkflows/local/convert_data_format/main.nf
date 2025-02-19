
include { UTILS_CONVERT_DATA_FORMAT_FOR_WORKFLOW  } from '../../modules/utils/convert_data_format_for_workflow'
workflow CONVERT_DATA_FORMAT {
    println("Input files:")
    params.convert_input =~ /\*/ ? file(params.convert_input).each{println "$it"} : println (file(params.convert_input))
    input_files = Channel.fromPath(params.convert_input)
    UTILS_CONVERT_DATA_FORMAT_FOR_WORKFLOW(input_files)

}
