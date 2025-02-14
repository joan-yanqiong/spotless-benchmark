nextflow.enable.dsl=2

// Deconvolution methods
include { MUSIC_RUN } from '../../../modules/local/music'
include { RCTD_RUN } from '../../../modules/local/rctd'
include { SPOTLIGHT_RUN } from '../../../modules/local/spotlight'
include { SPATIAL_DWLS_RUN } from '../../../modules/local/spatialdwls'
include { CELL2LOCATION_BUILD_MODEL} from '../../../modules/local/cell2location/build_model'
include {CELL2LOCATION_FIT_MODEL} from '../../../modules/local/cell2location/fit_model'

// include { STEREOSCOPE_BUILD_MODEL; STEREOSCOPE_FIT_MODEL } from '../../../modules/local/stereoscope'
// include { DESTVI_BUILD_MODEL; DESTVI_FIT_MODEL } from '../../../modules/local/destvi'
// include { runDSTG } from '../../../modules/local/dstg'
include { DSTG_CONVERT_SP_SEURAT_TO_MATRIX} from '../../../modules/local/dstg/convert_sp_seurat_to_matrix'
include { DSTG_CONVERT_REF_SEURAT_TO_MATRIX} from '../../../modules/local/dstg/convert_ref_seurat_to_matrix'

// include { NNLS_RUN } from '../../../modules/local/nnls'
// include { SEURAT_RUN } from '../../../modules/local/seurat'
// include { getCellComposition; TANGRAM_RUN } from '../../../modules/local/tangram'
// include { STRIDE_BUILD_MODEL; STRIDE_FIT_MODEL; STRIDE_RUN } from '../../../modules/local/stride'

// Helper functions
include { UTILS_CONVERT_BETWEEN_RDS_AND_H5AD as SC_OBJECT_CONVERT; UTILS_CONVERT_BETWEEN_RDS_AND_H5AD as SP_OBJECT_CONVERT } from '../../../modules/local/utils/convert_between_rds_and_h5ad'
include { UTILS_FORMAT_TSV_FILE as STEREOSCOPE_FORMAT; UTILS_FORMAT_TSV_FILE as CELL2LOCATION_FORMAT; UTILS_FORMAT_TSV_FILE as DESTVI_FORMAT; UTILS_FORMAT_TSV_FILE as DSTG_FORMAT;
          UTILS_FORMAT_TSV_FILE as TANGRAM_FORMAT; UTILS_FORMAT_TSV_FILE as STRIDE_FORMAT } from '../../../modules/local/utils/format_tsv_file'
include { UTILS_CREATE_DUMMY_FILE } from '../../../modules/local/utils/create_dummy_file'

workflow RUN_DECONVOLUTION {
    take:
        ch_ref_obj
        ch_sp_obj
        // sc_input_type
        // sp_input_type

    main:
        // String matching to check which method to run
        all_methods = "music,rctd,spatialdwls,spotlight,stereoscope,cell2location,destvi,dstg,nnls,seurat,tangram,stride"
        r_methods = ["music", "rctd", "spatialdwls", "spotlight", "dstg", "nnls", "seurat"]
        python_methods = ["stereoscope","cell2location","destvi","tangram","stride"]

        methods = ( params.methods.toLowerCase() ==~ /all/ ? all_methods : params.methods.toLowerCase() )
        methods_list = Arrays.asList(methods.split(','))
        
        output_ch = Channel.empty() // collect output channels

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

        ch_sp_obj_only = ch_sp_obj.map {
            it -> it[1]
        }

        ch_sp_obj
            .map {
            // it[0] sample_id, it[1] filepath
            // use it[1] to create key for merging
            sample_id, filepath -> [filepath.simpleName, sample_id, filepath]
            }
            // Out: [0] key, [1] sample_id, [2] filepath (rds/h5ad), [3] filepath (rds/h5ad)
            .join(
                SP_OBJECT_CONVERT(ch_sp_obj_only)
                    .map {
                        it -> [it.simpleName, it]
                    }     
            )
            .multiMap {
                _key, sample_id, filepath1, filepath2 -> 
                h5ad: [sample_id, (filepath1.extension == "h5ad" ? filepath1 : filepath2)]
                rds: [sample_id, (filepath1.extension == "rds" ? filepath1 : filepath2)]
            }
            .set { ch_sp_paired_obj }   

        
        ch_ref_paired_obj.rds.combine(ch_sp_paired_obj.rds)
            .multiMap { ref_obj, sample_id, sp_obj -> 
                    ref_obj: ref_obj
                    sp_input: [sample_id, sp_obj]
            }
            .set { ch_r_methods_input }

        // // R methods
        if ( !methods_list.disjoint(r_methods) ) {

            // if ( methods =~ /music/ ){
            //     MUSIC_RUN(ch_r_methods_input.ref_obj, ch_r_methods_input.sp_input)
            //     output_ch = output_ch.mix(MUSIC_RUN.out.txt)
            // }
            if ( methods =~ /rctd/ ){
                RCTD_RUN(ch_r_methods_input.ref_obj, ch_r_methods_input.sp_input)
                output_ch = output_ch.mix(RCTD_RUN.out.txt)
            }
            
            if ( methods =~ /spotlight/ ){
                SPOTLIGHT_RUN(ch_r_methods_input.ref_obj, ch_r_methods_input.sp_input)
                output_ch = output_ch.mix(SPOTLIGHT_RUN.out)
            }

            if ( methods =~ /spatialdwls/ ){
                SPATIAL_DWLS_RUN(ch_r_methods_input.ref_obj, ch_r_methods_input.sp_input)
                output_ch = output_ch.mix(SPATIAL_DWLS_RUN.out.txt)
            }

            // if ( methods =~ /dstg/ ){
            //     DSTG_CONVERT_REF_SEURAT_TO_MATRIX(ch_ref_paired_obj.rds)
            //     DSTG_CONVERT_SP_SEURAT_TO_MATRIX(ch_sp_paired_obj.rds)

            //     // DSTG_CONVERT_REF_SEURAT_TO_MATRIX.out.rds.combine(
            //     //     DSTG_CONVERT_SP_SEURAT_TO_MATRIX.out.rds
            //     // )
            //     // .multiMap(
                    
            //     //     ref_input: 
            //     // )

            //     // runDSTG(pair_input_ch)
            //     DSTG_FORMAT(runDSTG.out.txt) 
            //     // output_ch = output_ch.mix(DSTG_FORMAT.out)
            // }

        //     // if ( methods =~ /nnls/ ){
        //     //     NNLS_RUN(pair_input_ch)
        //     //     output_ch = output_ch.mix(NNLS_RUN.out)
        //     // }

        //     // if ( methods =~ /seurat/ ){
        //     //     SEURAT_RUN(pair_input_ch)
        //     //     output_ch = output_ch.mix(SEURAT_RUN.out)
        //     // }
        }
        // Python methods
        // First check if there are python methods in the input params
        // before performing conversion of data to h5ad
        
        if ( !methods_list.disjoint(python_methods) ){

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

            if ( methods =~ /cell2location/ ) {
                CELL2LOCATION_BUILD_MODEL(ch_ref_paired_obj.h5ad)

                // Repeat model output for each spatial file
                CELL2LOCATION_BUILD_MODEL.out.h5ad.combine(ch_sp_paired_obj.h5ad)
                .multiMap { model_sc_file, sample_id, sp_file_h5ad ->
                            model: model_sc_file
                            sp_input: [sample_id, sp_file_h5ad] }
                .set{ ch_c2l_input }

                CELL2LOCATION_FIT_MODEL(
                    ch_c2l_input.model,
                    ch_c2l_input.sp_input
                )
                CELL2LOCATION_FORMAT(CELL2LOCATION_FIT_MODEL.out)
                output_ch = output_ch.mix(CELL2LOCATION_FORMAT.out)
            }

            // if ( methods =~ /destvi/ ) {
                
            //     DESTVI_BUILD_MODEL(sc_input_conv, sp_input_pair.first())
                
            //     // Repeat model output for each spatial file
            //     DESTVI_BUILD_MODEL.out.combine(sp_input_pair)
            //     .multiMap { h5ad_file, model_file, sp_file_h5ad, sp_file_rds ->
            //                 model: tuple h5ad_file, model_file
            //                 sp_input: tuple sp_file_h5ad, sp_file_rds }
            //     .set{ destvi_combined_ch }

            //     DESTVI_FIT_MODEL(destvi_combined_ch.sp_input,
            //                    destvi_combined_ch.model)
            //     DESTVI_FORMAT(DESTVI_FIT_MODEL.out) 
            //     output_ch = output_ch.mix(DESTVI_FORMAT.out)
            // }

            // if ( methods =~ /tangram/ ){
            //     // For constrained mode, need to run extra function
            //     if (params.deconv_args.tangram =~ /-m[ ]*constrained/) {
            //         // Get cell counts from spatial file (or segment)
            //         getCellComposition(sp_input_pair)

            //         // Repeat cell composition and spatial file for each single-cell file
            //         getCellComposition.out.combine(sc_input_conv)
            //         .multiMap { cell_count_file, sp_file_h5ad, sp_file_rds, sc_file ->
            //                     sc_input: sc_file
            //                     sp_input: tuple sp_file_h5ad, sp_file_rds
            //                     cell_count: cell_count_file }
            //         .set{ tangram_combined_ch }
                    
            //         TANGRAM_RUN(tangram_combined_ch.sc_input,  tangram_combined_ch.sp_input,
            //                    tangram_combined_ch.cell_count)
            //     } else {
            //         // Only take single-cell and spatial file
            //         sc_input_conv.combine(sp_input_pair)
            //         .multiMap { sc_file, sp_file_h5ad, sp_file_rds ->
            //                     sc_input: sc_file
            //                     sp_input: tuple sp_file_h5ad, sp_file_rds}
            //         .set{ tangram_combined_ch }

            //         // Cell count file as empty list
            //         TANGRAM_RUN(tangram_combined_ch.sc_input,  tangram_combined_ch.sp_input, [])

            //     }

            //     TANGRAM_FORMAT(TANGRAM_RUN.out)
            //     output_ch = output_ch.mix(TANGRAM_FORMAT.out)

            // }

            // if ( methods =~ /stride/ ) {
                
            //     /*
            //     STRIDE_BUILD_MODEL(sc_input_conv)

            //     // Repeat model output for each spatial file
            //     STRIDE_BUILD_MODEL.out.combine(sp_input_pair)
            //     .multiMap { zip_file, txt_files, sp_file_h5ad, sp_file_rds ->
            //                 model: tuple zip_file, txt_files
            //                 sp_input: tuple sp_file_h5ad, sp_file_rds }
            //     .set{ stride_combined_ch }
                
            //     STRIDE_FIT_MODEL(stride_combined_ch.sp_input,
            //                    stride_combined_ch.model)
            //     STRIDE_FORMAT(STRIDE_FIT_MODEL.out) 
            //     */
                
                
            //     STRIDE_RUN(sc_input_conv.combine(sp_input_pair))
            //     STRIDE_FORMAT(STRIDE_RUN.out)
                
            //     output_ch = output_ch.mix(STRIDE_FORMAT.out)
            // }
        }

    emit:
        output_ch
}

workflow {

    if (!(params.mode ==~ /run_dataset/)){
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
        params.sp_input =~ /\*/ ? file(params.sp_input).each{println "$it"} : println (file(params.sp_input))
    }

    ch_ref_obj = Channel.fromPath(params.sc_input) // Can only have 1 file
    ch_sp_obj = Channel.fromPath(params.sp_input) // Can have one or more files

    RUN_DECONVOLUTION(ch_ref_obj, ch_sp_obj)
}
