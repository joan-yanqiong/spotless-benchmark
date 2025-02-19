nextflow.enable.dsl = 2

// Deconvolution methods
// include { MUSIC_RUN } from '../../../modules/local/music'
// include { MUSIC_DOWNSAMPLE_SCRNASEQ_REF } from '../../../modules/local/music/downsample_scrnaseq_ref'
// include { MUSIC_CHECK_SCRNASEQ_REF_SIZE } from '../../../modules/local/music/check_scrnaseq_ref_size'
include { RCTD_RUN                                                } from '../../../modules/local/rctd'
include { SPOTLIGHT_RUN                                           } from '../../../modules/local/spotlight'
include { SPATIALDWLS_RUN                                         } from '../../../modules/local/spatialdwls'
include { CELL2LOCATION_BUILD_MODEL                               } from '../../../modules/local/cell2location/build_model'
include { CELL2LOCATION_FIT_MODEL                                 } from '../../../modules/local/cell2location/fit_model'

// include { STEREOSCOPE_BUILD_MODEL; STEREOSCOPE_FIT_MODEL } from '../../../modules/local/stereoscope'
// include { DESTVI_BUILD_MODEL; DESTVI_FIT_MODEL } from '../../../modules/local/destvi'
// include { runDSTG } from '../../../modules/local/dstg'
// include { DSTG_CONVERT_SP_SEURAT_TO_MATRIX} from '../../../modules/local/dstg/convert_sp_seurat_to_matrix'
// include { DSTG_CONVERT_REF_SEURAT_TO_MATRIX} from '../../../modules/local/dstg/convert_ref_seurat_to_matrix'

// include { NNLS_RUN } from '../../../modules/local/nnls'
// include { SEURAT_RUN } from '../../../modules/local/seurat'
// include { getCellComposition; TANGRAM_RUN } from '../../../modules/local/tangram'
// include { STRIDE_BUILD_MODEL; STRIDE_FIT_MODEL; STRIDE_RUN } from '../../../modules/local/stride'

// Helper functions
include { UTILS_CONVERT_BETWEEN_RDS_AND_H5AD as SC_OBJECT_CONVERT ; UTILS_CONVERT_BETWEEN_RDS_AND_H5AD as SP_OBJECT_CONVERT } from '../../../modules/local/utils/convert_between_rds_and_h5ad'
include {
    UTILS_FORMAT_TSV_FILE as STEREOSCOPE_FORMAT ;
    UTILS_FORMAT_TSV_FILE as CELL2LOCATION_FORMAT ;
    UTILS_FORMAT_TSV_FILE as DESTVI_FORMAT ;
    UTILS_FORMAT_TSV_FILE as DSTG_FORMAT ;
    UTILS_FORMAT_TSV_FILE as TANGRAM_FORMAT ;
    UTILS_FORMAT_TSV_FILE as STRIDE_FORMAT
} from '../../../modules/local/utils/format_tsv_file'
include { UTILS_CREATE_DUMMY_FILE                                 } from '../../../modules/local/utils/create_dummy_file'

workflow {

    if (!(params.mode ==~ /run_dataset/)) {
        throw new Exception("Error: can only run this with the 'run_dataset' mode")
    }

    // RUN ON YOUR OWN DATA
    println("Running the pipeline on the provided data...")

    // Print inputs (the timing isn't right with with view(), so do this instead)
    // Although view() has the advantage that it gives the absolute path (ch_ref_obj.view())
    if (params.verbose) {
        println("Single-cell reference dataset:")
        println(file(params.sc_input))

        println("\nSpatial dataset(s):")
        // With glob pattern, there will be multiple files
        params.sp_input =~ /\*/ ? file(params.sp_input).each { println("${it}") } : println(file(params.sp_input))
    }

    ch_ref_obj = Channel.fromPath(params.sc_input)
    // Can only have 1 file
    ch_sp_obj = Channel.fromPath(params.sp_input)
    // Can have one or more files

    RUN_DECONVOLUTION(ch_ref_obj, ch_sp_obj)
}

workflow RUN_DECONVOLUTION {
    take:
    ch_ref_obj
    ch_sp_obj

    main:
    // String matching to check which method to run
    all_methods = "music,rctd,spatialdwls,spotlight,stereoscope,cell2location,destvi,dstg,nnls,seurat,tangram,stride"
    r_methods = ["music", "rctd", "spatialdwls", "spotlight", "dstg", "nnls", "seurat"]
    python_methods = ["stereoscope", "cell2location", "destvi", "tangram", "stride"]

    methods = (params.methods.toLowerCase() ==~ /all/ ? all_methods : params.methods.toLowerCase())
    methods_list = Arrays.asList(methods.split(','))

    output_ch = Channel.empty()
    // collect output channels

    // ---- Prepare scRNAseq reference ----
    SC_OBJECT_CONVERT(ch_ref_obj)
        .combine(ch_ref_obj)
        .flatten()
        .branch { it ->
            rds: it.extension == "rds"
            return it
            h5ad: it.extension == "h5ad"
            return it
            other: true
        }
        .set { ch_ref_paired_obj }

    ch_sp_obj_only = ch_sp_obj.map { it ->
        it[1]
    }

    ch_sp_obj
        .map { meta, filepath ->
            [filepath.simpleName, meta, filepath]
        }
        .join(
            SP_OBJECT_CONVERT(ch_sp_obj_only).map { it ->
                [it.simpleName, it]
            }
        )
        .multiMap { _key, meta, filepath1, filepath2 ->
            h5ad: [meta, (filepath1.extension == "h5ad" ? filepath1 : filepath2)]
            rds: [meta, (filepath1.extension == "rds" ? filepath1 : filepath2)]
        }
        .set { ch_sp_paired_obj }


    ch_ref_paired_obj.rds
        .combine(ch_sp_paired_obj.rds)
        .multiMap { ref_obj, meta, sp_obj ->
            ref_obj: ref_obj
            sp_input: [meta, sp_obj]
        }
        .set { ch_r_methods_input }

    // // R methods
    if (!methods_list.disjoint(r_methods)) {

        // if ( methods =~ /music/ ){
        //     MUSIC_CHECK_SCRNASEQ_REF_SIZE(ch_ref_paired_obj.rds)

        //     MUSIC_CHECK_SCRNASEQ_REF_SIZE.out.rds.ifEmpty(
        //         ch_ref_paired_obj.rds).view()
        // }

        // MUSIC_RUN(ch_r_methods_input.ref_obj, ch_r_methods_input.sp_input)
        // output_ch = output_ch.mix(MUSIC_RUN.out.txt)

        if (methods =~ /rctd/) {
            RCTD_RUN(ch_r_methods_input.ref_obj, ch_r_methods_input.sp_input)
            output_ch = output_ch.mix(RCTD_RUN.out.txt)
        }

        if (methods =~ /spotlight/) {
            SPOTLIGHT_RUN(ch_r_methods_input.ref_obj, ch_r_methods_input.sp_input)
            ch_r_methods_input.ref_obj.view()
            output_ch = output_ch.mix(SPOTLIGHT_RUN.out)
        }

        if (methods =~ /spatialdwls/) {
            SPATIALDWLS_RUN(ch_r_methods_input.ref_obj, ch_r_methods_input.sp_input)
            output_ch = output_ch.mix(SPATIALDWLS_RUN.out.txt)
        }
    }
    // Python methods
    // First check if there are python methods in the input params
    // before performing conversion of data to h5ad

    if (!methods_list.disjoint(python_methods)) {

        // if ( methods =~ /stereoscope/ ) {
        //     STEREOSCOPE_BUILD_MODEL(sc_input_conv)

        //     // Repeat model output for each spatial file
        //     STEREOSCOPE_BUILD_MODEL.out.combine(sp_input_pair)
        //     .multiMap { r_file, logits_file, sp_file_h5ad, sp_file_rds ->
        //                 model: tuple r_file, logits_file
        //                 sp_input: tuple sp_file_h5ad, sp_file_rds }
        //     .set{ stereo_combined_ch }

        //     STEREOSCOPE_FIT_MODEL(stereo_combined_ch.sp_input,
        //                         stereo_combined_ch.model)
        //     STEREOSCOPE_FORMAT(STEREOSCOPE_FIT_MODEL.out) 
        //     output_ch = output_ch.mix(STEREOSCOPE_FORMAT.out)
        // }

        if (methods =~ /cell2location/) {
            CELL2LOCATION_BUILD_MODEL(ch_ref_paired_obj.h5ad)

            // Repeat model output for each spatial file
            CELL2LOCATION_BUILD_MODEL.out.h5ad
                .combine(ch_sp_paired_obj.h5ad)
                .multiMap { model_sc_file, meta, sp_file_h5ad ->
                    model: model_sc_file
                    sp_input: [meta, sp_file_h5ad]
                }
                .set { ch_c2l_input }

            CELL2LOCATION_FIT_MODEL(
                ch_c2l_input.model,
                ch_c2l_input.sp_input,
            )
            CELL2LOCATION_FORMAT(CELL2LOCATION_FIT_MODEL.out)
            output_ch = output_ch.mix(CELL2LOCATION_FORMAT.out)
        }
    }

    emit:
    output_ch
}
