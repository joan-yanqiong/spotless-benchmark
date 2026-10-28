include { COMPUTE_METRICS } from '../../../modules/local/compute_metrics'

// In case metrics.R is updated, just rerun metrics computation
// Give sp_input and methods as parameters
workflow {
    sp_input_ch = Channel.fromPath(params.sp_input)
    
    all_methods = "music,rctd,spatialdwls,spotlight,stereoscope,cell2location,destvi,dstg,nnls,seurat,tangram,stride"
    methods = ( params.methods ==~ /all/ ? all_methods : params.methods )
    methods_list = Channel.from(Arrays.asList(methods.split(',')))

    // Get directory of the proportions file
    input_ch = methods_list.combine(sp_input_ch)
                .multiMap { method, sp_input ->
                    inputs: tuple (method, \
                    "${params.outdir.props}/${file(sp_input).getSimpleName().replaceFirst(/_[a-z]{3}[0-9]+/, "")}/proportions_${method}_${file(sp_input).getSimpleName()}${params.runID_props}",\
                    sp_input)
                }
                
    COMPUTE_METRICS(input_ch.inputs)

}