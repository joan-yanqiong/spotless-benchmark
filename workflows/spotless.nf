/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_spotless_pipeline'

include { RUN_DECONVOLUTION      } from '../subworkflows/local/deconvolution/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPOTLESS {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_ref_obj     // value channel: [ path(sc_input) ]

    main:

    ch_versions = Channel.empty()

    // // Define constants
    // silver_standards = ((1..7)*9).sort().withIndex().collect{ it, i -> "silver_standard_$it-${ -> i%9+1}".toString() }
    // silver_standards_coarse = (1..7).collect{ it -> "silver_standard_$it".toString()}
    // gold_standards = (1..3).collect{ it -> "gold_standard_$it".toString()}
    // all_standards = gold_standards + silver_standards + silver_standards_coarse + ["silver_standard_1-10", "silver_standard_1-11"]
    all_modes = ["run_standard", "run_dataset", "generate_and_run"]

    if (!all_modes.contains(params.mode)) {
        throw new Exception("Error: please enter --mode as 'run_standard', 'run_dataset', or 'generate_and_run'")
    }

    // // RUN STANDARD PIPELINE
    // if (params.mode ==~ /run_standard/ ){
    //     // Check pipeline name
    //     if (!all_standards.contains(params.standard)){
    //         throw new Exception("Error: standard '$params.standard' does not exist.")
    //     }
    //     // Check if it was run with config file
    //     if (params.sc_input =~ /test_sc_data/ && params.sp_input =~ /test_sp_data/){
    //         throw new Exception("Error: default input detected, please run with '-c data/pipeline.config'")
    //     }
    //     println("Running the pipeline on $params.standard...")
    // } 
    // // RUN ON YOUR OWN DATA
    // else if (params.mode ==~ /run_dataset/ ) {
    //     println("Running the pipeline on the provided data...")
    // }
    // // GENERATE SYNTHETIC DATA AND RUN
    // else if (params.mode ==~ /generate_and_run/) {
    //     println("Generating synthetic data with synthspot...")
    //     GENERATE_SYNTHETIC_DATA()
    //     println("Synthetic data will be copied to ${params.outdir.synthspot}")
    // }

    // // Print inputs (the timing isn't right with with view(), so do this instead)
    // // Although view() has the advantage that it gives the absolute path (sc_input_ch.view())
    // if (!(params.mode ==~ /generate_and_run/)) {
    //     if (params.verbose) {
    //         println("Single-cell reference dataset:")
    //         println(params.mode ==~ /run_standard/ ? file(params.sc_input)[0]: file(params.sc_input))

    //         println("\nSpatial dataset(s):") 
    //         // With glob pattern, there will be multiple files
    //         params.sp_input =~ /\*/ ? file(params.sp_input).each{println "$it"} : println (file(params.sp_input))
    //     }
    // }

    // sc_input_type = params.sc_input.split("\\.").last() =~ /h5/ ? "h5ad" : "rds"
    // sp_input_type = ch_samplesheet[1].split("\\.").last() =~ /h5/ ? "h5ad" : "rds"

    // println("The single-cell data is of ${sc_input_type} format and the spatial data is of ${sp_input_type} format.")

    // // Print arguments
    // println("\nParameters:")
    // params.deconv_args.each{ method, argu  -> println "${method}: ${argu}" }

    // // if (sp_input_type ==~ /h5ad/ && !params.skip_metrics){
    // //     throw new Exception("Metric calculation is not supported for spatial h5ad files as input, please use --skip_metrics if you wish to run the pipeline.")
    // // }


    // // sc_input_ch
    // // // // Can have one or more files
    // // sp_input_ch = (params.mode ==~ /generate_and_run/ ? GENERATE_SYNTHETIC_DATA.out : ch_samplesheet)

    RUN_DECONVOLUTION(ch_ref_obj, ch_samplesheet)

    // if (!params.skip_metrics){
    //     COMPUTE_METRICS(RUN_DECONVOLUTION.out)
    // }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'spotless_software_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
}
