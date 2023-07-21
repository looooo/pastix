/**
 *
 * @file check_options.c
 *
 * This file is generated automatically. If you want to modify it, modify
 * ${PASTIX_HOME}/tools/gen_param/pastix_[iparm/dparm/enums].py and run
 * ${PASTIX_HOME}/tools/gen_param/gen_parm_files.py ${PASTIX_HOME}.
 *
 * @copyright 2004-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Esragul Korkmaz
 * @author Gregoire Pichon
 * @author Tony Delarue
 * @date 2023-04-19
 *
 */
#include "common.h"

static inline int
pastix_verbose_check_value( pastix_verbose_t value )
{
    if( (value == PastixVerboseNot) ||
        (value == PastixVerboseNo) ||
        (value == PastixVerboseYes) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_io_check_value( pastix_io_t value )
{
    if( (value == PastixIONo) ||
        (value == PastixIOLoad) ||
        (value == PastixIOSave) ||
        (value == PastixIOLoadGraph) ||
        (value == PastixIOSaveGraph) ||
        (value == PastixIOLoadCSC) ||
        (value == PastixIOSaveCSC) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_trace_check_value( pastix_trace_t value )
{
    if( (value == PastixTraceNumfact) ||
        (value == PastixTraceSolve) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_ordering_check_value( pastix_ordering_t value )
{
    if( (value == PastixOrderScotch) ||
        (value == PastixOrderMetis) ||
        (value == PastixOrderPersonal) ||
        (value == PastixOrderPtScotch) ||
        (value == PastixOrderParMetis) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_split_check_value( pastix_split_t value )
{
    if( (value == PastixSplitNot) ||
        (value == PastixSplitKway) ||
        (value == PastixSplitKwayProjections) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_factotype_check_value( pastix_factotype_t value )
{
    if( (value == PastixFactLLH) ||
        (value == PastixFactLDLT) ||
        (value == PastixFactLU) ||
        (value == PastixFactLLT) ||
        (value == PastixFactLDLH) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_fact_mode_check_value( pastix_fact_mode_t value )
{
    if( (value == PastixFactModeLocal) ||
        (value == PastixFactModeSchur) ||
        (value == PastixFactModeBoth) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_trans_check_value( pastix_trans_t value )
{
    if( (value == PastixNoTrans) ||
        (value == PastixTrans) ||
        (value == PastixConjTrans) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_solv_mode_check_value( pastix_solv_mode_t value )
{
    if( (value == PastixSolvModeLocal) ||
        (value == PastixSolvModeInterface) ||
        (value == PastixSolvModeSchur) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_refine_check_value( pastix_refine_t value )
{
    if( (value == PastixRefineGMRES) ||
        (value == PastixRefineCG) ||
        (value == PastixRefineSR) ||
        (value == PastixRefineBiCGSTAB) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_scheduler_check_value( pastix_scheduler_t value )
{
    if( (value == PastixSchedSequential) ||
        (value == PastixSchedStatic) ||
        (value == PastixSchedParsec) ||
        (value == PastixSchedStarPU) ||
        (value == PastixSchedDynamic) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_compress_when_check_value( pastix_compress_when_t value )
{
    if( (value == PastixCompressNever) ||
        (value == PastixCompressWhenBegin) ||
        (value == PastixCompressWhenEnd) ||
        (value == PastixCompressWhenDuring) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_compress_method_check_value( pastix_compress_method_t value )
{
    if( (value == PastixCompressMethodSVD) ||
        (value == PastixCompressMethodPQRCP) ||
        (value == PastixCompressMethodRQRCP) ||
        (value == PastixCompressMethodTQRCP) ||
        (value == PastixCompressMethodRQRRT) ||
        (value == PastixCompressMethodNbr) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_compress_ortho_check_value( pastix_compress_ortho_t value )
{
    if( (value == PastixCompressOrthoCGS) ||
        (value == PastixCompressOrthoQR) ||
        (value == PastixCompressOrthoPartialQR) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_task_check_value( pastix_task_t value )
{
    if( (value == PastixTaskInit) ||
        (value == PastixTaskOrdering) ||
        (value == PastixTaskSymbfact) ||
        (value == PastixTaskAnalyze) ||
        (value == PastixTaskNumfact) ||
        (value == PastixTaskSolve) ||
        (value == PastixTaskRefine) ||
        (value == PastixTaskClean) ) {
        return 0;
    }
    return 1;
}

static inline int
pastix_coeftype_check_value( pastix_coeftype_t value )
{
    if( (value == PastixPattern) ||
        (value == PastixFloat) ||
        (value == PastixDouble) ||
        (value == PastixComplex32) ||
        (value == PastixComplex64) ) {
        return 0;
    }
    return 1;
}

static inline int
iparm_verbose_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_verbose_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_VERBOSE: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_io_strategy_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_io_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_IO_STRATEGY: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_produce_stats_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_PRODUCE_STATS] */
    (void)iparm;
    return 0;
}

static inline int
iparm_trace_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_trace_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_TRACE: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_mc64_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_MC64] */
    (void)iparm;
    return 0;
}

static inline int
iparm_ordering_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_ordering_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_ORDERING: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_ordering_default_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_ORDERING_DEFAULT] */
    (void)iparm;
    return 0;
}

static inline int
iparm_scotch_mt_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_SCOTCH_MT] */
    (void)iparm;
    return 0;
}

static inline int
iparm_scotch_switch_level_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_SCOTCH_SWITCH_LEVEL] */
    (void)iparm;
    return 0;
}

static inline int
iparm_scotch_cmin_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_SCOTCH_CMIN] */
    (void)iparm;
    return 0;
}

static inline int
iparm_scotch_cmax_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_SCOTCH_CMAX] */
    (void)iparm;
    return 0;
}

static inline int
iparm_scotch_frat_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_SCOTCH_FRAT] */
    (void)iparm;
    return 0;
}

static inline int
iparm_metis_ctype_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_METIS_CTYPE] */
    (void)iparm;
    return 0;
}

static inline int
iparm_metis_rtype_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_METIS_RTYPE] */
    (void)iparm;
    return 0;
}

static inline int
iparm_metis_no2hop_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_METIS_NO2HOP] */
    (void)iparm;
    return 0;
}

static inline int
iparm_metis_nseps_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_METIS_NSEPS] */
    (void)iparm;
    return 0;
}

static inline int
iparm_metis_niter_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_METIS_NITER] */
    (void)iparm;
    return 0;
}

static inline int
iparm_metis_ufactor_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_METIS_UFACTOR] */
    (void)iparm;
    return 0;
}

static inline int
iparm_metis_compress_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_METIS_COMPRESS] */
    (void)iparm;
    return 0;
}

static inline int
iparm_metis_ccorder_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_METIS_CCORDER] */
    (void)iparm;
    return 0;
}

static inline int
iparm_metis_pfactor_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_METIS_PFACTOR] */
    (void)iparm;
    return 0;
}

static inline int
iparm_metis_seed_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_METIS_SEED] */
    (void)iparm;
    return 0;
}

static inline int
iparm_metis_dbglvl_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_METIS_DBGLVL] */
    (void)iparm;
    return 0;
}

static inline int
iparm_amalgamation_lvlblas_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_AMALGAMATION_LVLBLAS] */
    (void)iparm;
    return 0;
}

static inline int
iparm_amalgamation_lvlcblk_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_AMALGAMATION_LVLCBLK] */
    (void)iparm;
    return 0;
}

static inline int
iparm_reordering_split_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_REORDERING_SPLIT] */
    (void)iparm;
    return 0;
}

static inline int
iparm_reordering_stop_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_REORDERING_STOP] */
    (void)iparm;
    return 0;
}

static inline int
iparm_splitting_strategy_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_split_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_SPLITTING_STRATEGY: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_splitting_levels_projections_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_SPLITTING_LEVELS_PROJECTIONS] */
    (void)iparm;
    return 0;
}

static inline int
iparm_splitting_levels_kway_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_SPLITTING_LEVELS_KWAY] */
    (void)iparm;
    return 0;
}

static inline int
iparm_splitting_projections_depth_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_SPLITTING_PROJECTIONS_DEPTH] */
    (void)iparm;
    return 0;
}

static inline int
iparm_splitting_projections_distance_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_SPLITTING_PROJECTIONS_DISTANCE] */
    (void)iparm;
    return 0;
}

static inline int
iparm_splitting_projections_width_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_SPLITTING_PROJECTIONS_WIDTH] */
    (void)iparm;
    return 0;
}

static inline int
iparm_min_blocksize_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_MIN_BLOCKSIZE] */
    (void)iparm;
    return 0;
}

static inline int
iparm_max_blocksize_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_MAX_BLOCKSIZE] */
    (void)iparm;
    return 0;
}

static inline int
iparm_tasks2d_level_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_TASKS2D_LEVEL] */
    (void)iparm;
    return 0;
}

static inline int
iparm_tasks2d_width_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_TASKS2D_WIDTH] */
    (void)iparm;
    return 0;
}

static inline int
iparm_allcand_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_ALLCAND] */
    (void)iparm;
    return 0;
}

static inline int
iparm_incomplete_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_INCOMPLETE] */
    (void)iparm;
    return 0;
}

static inline int
iparm_level_of_fill_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_LEVEL_OF_FILL] */
    (void)iparm;
    return 0;
}

static inline int
iparm_factorization_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_factotype_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_FACTORIZATION: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_facto_look_side_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_FACTO_LOOK_SIDE] */
    (void)iparm;
    return 0;
}

static inline int
iparm_free_cscuser_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_FREE_CSCUSER] */
    (void)iparm;
    return 0;
}

static inline int
iparm_schur_fact_mode_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_fact_mode_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_SCHUR_FACT_MODE: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_transpose_solve_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_trans_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_TRANSPOSE_SOLVE: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_schur_solv_mode_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_solv_mode_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_SCHUR_SOLV_MODE: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_applyperm_ws_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_APPLYPERM_WS] */
    (void)iparm;
    return 0;
}

static inline int
iparm_refinement_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_refine_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_REFINEMENT: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_itermax_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_ITERMAX] */
    (void)iparm;
    return 0;
}

static inline int
iparm_gmres_im_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_GMRES_IM] */
    (void)iparm;
    return 0;
}

static inline int
iparm_scheduler_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_scheduler_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_SCHEDULER: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_thread_nbr_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_THREAD_NBR] */
    (void)iparm;
    return 0;
}

static inline int
iparm_autosplit_comm_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_AUTOSPLIT_COMM] */
    (void)iparm;
    return 0;
}

static inline int
iparm_gpu_nbr_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_GPU_NBR] */
    (void)iparm;
    return 0;
}

static inline int
iparm_gpu_memory_percentage_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_GPU_MEMORY_PERCENTAGE] */
    (void)iparm;
    return 0;
}

static inline int
iparm_gpu_memory_block_size_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_GPU_MEMORY_BLOCK_SIZE] */
    (void)iparm;
    return 0;
}

static inline int
iparm_global_allocation_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_GLOBAL_ALLOCATION] */
    (void)iparm;
    return 0;
}

static inline int
iparm_compress_min_width_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_COMPRESS_MIN_WIDTH] */
    (void)iparm;
    return 0;
}

static inline int
iparm_compress_min_height_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_COMPRESS_MIN_HEIGHT] */
    (void)iparm;
    return 0;
}

static inline int
iparm_compress_when_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_compress_when_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_COMPRESS_WHEN: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_compress_method_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_compress_method_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_COMPRESS_METHOD: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_compress_ortho_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_compress_ortho_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_COMPRESS_ORTHO: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_compress_reltol_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_COMPRESS_RELTOL] */
    (void)iparm;
    return 0;
}

static inline int
iparm_compress_preselect_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_COMPRESS_PRESELECT] */
    (void)iparm;
    return 0;
}

static inline int
iparm_compress_iluk_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_COMPRESS_ILUK] */
    (void)iparm;
    return 0;
}

static inline int
iparm_mixed_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_MIXED] */
    (void)iparm;
    return 0;
}

static inline int
iparm_ftz_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_FTZ] */
    (void)iparm;
    return 0;
}

static inline int
iparm_modify_parameter_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_MODIFY_PARAMETER] */
    (void)iparm;
    return 0;
}

static inline int
iparm_start_task_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_task_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_START_TASK: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_float_check_value( pastix_int_t iparm )
{
    int rc;
    rc = pastix_coeftype_check_value( iparm );
    if ( rc == 1 ) {
        fprintf(stderr, "IPARM_FLOAT: The value is incorrect\n");
    }
    return rc;
}

static inline int
iparm_mtx_type_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_MTX_TYPE] */
    (void)iparm;
    return 0;
}

static inline int
iparm_dof_nbr_check_value( pastix_int_t iparm )
{
    /* TODO : Check range iparm[IPARM_DOF_NBR] */
    (void)iparm;
    return 0;
}

static inline int
dparm_epsilon_refinement_check_value( double dparm )
{
    /* TODO : Check range dparm[DPARM_EPSILON_REFINEMENT] */
    (void)dparm;
    return 0;
}

static inline int
dparm_epsilon_magn_ctrl_check_value( double dparm )
{
    /* TODO : Check range dparm[DPARM_EPSILON_MAGN_CTRL] */
    (void)dparm;
    return 0;
}

static inline int
dparm_compress_tolerance_check_value( double dparm )
{
    /* TODO : Check range dparm[DPARM_COMPRESS_TOLERANCE] */
    (void)dparm;
    return 0;
}

static inline int
dparm_compress_min_ratio_check_value( double dparm )
{
    /* TODO : Check range dparm[DPARM_COMPRESS_MIN_RATIO] */
    (void)dparm;
    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Check the values of the iparm.
 *
 *******************************************************************************
 *
 * @param[in] iparm
 *          The iparm options array.
 *
 *******************************************************************************
 *
 * @return The amount of incorrect values in the iparm array.
 *
 *******************************************************************************/
int
iparm_check_values( const pastix_int_t *iparm )
{
    int error = 0;
    error += iparm_verbose_check_value( iparm[IPARM_VERBOSE] );
    error += iparm_io_strategy_check_value( iparm[IPARM_IO_STRATEGY] );
    error += iparm_produce_stats_check_value( iparm[IPARM_PRODUCE_STATS] );
    error += iparm_trace_check_value( iparm[IPARM_TRACE] );
    error += iparm_mc64_check_value( iparm[IPARM_MC64] );
    error += iparm_ordering_check_value( iparm[IPARM_ORDERING] );
    error += iparm_ordering_default_check_value( iparm[IPARM_ORDERING_DEFAULT] );
    error += iparm_scotch_mt_check_value( iparm[IPARM_SCOTCH_MT] );
    error += iparm_scotch_switch_level_check_value( iparm[IPARM_SCOTCH_SWITCH_LEVEL] );
    error += iparm_scotch_cmin_check_value( iparm[IPARM_SCOTCH_CMIN] );
    error += iparm_scotch_cmax_check_value( iparm[IPARM_SCOTCH_CMAX] );
    error += iparm_scotch_frat_check_value( iparm[IPARM_SCOTCH_FRAT] );
    error += iparm_metis_ctype_check_value( iparm[IPARM_METIS_CTYPE] );
    error += iparm_metis_rtype_check_value( iparm[IPARM_METIS_RTYPE] );
    error += iparm_metis_no2hop_check_value( iparm[IPARM_METIS_NO2HOP] );
    error += iparm_metis_nseps_check_value( iparm[IPARM_METIS_NSEPS] );
    error += iparm_metis_niter_check_value( iparm[IPARM_METIS_NITER] );
    error += iparm_metis_ufactor_check_value( iparm[IPARM_METIS_UFACTOR] );
    error += iparm_metis_compress_check_value( iparm[IPARM_METIS_COMPRESS] );
    error += iparm_metis_ccorder_check_value( iparm[IPARM_METIS_CCORDER] );
    error += iparm_metis_pfactor_check_value( iparm[IPARM_METIS_PFACTOR] );
    error += iparm_metis_seed_check_value( iparm[IPARM_METIS_SEED] );
    error += iparm_metis_dbglvl_check_value( iparm[IPARM_METIS_DBGLVL] );
    error += iparm_amalgamation_lvlblas_check_value( iparm[IPARM_AMALGAMATION_LVLBLAS] );
    error += iparm_amalgamation_lvlcblk_check_value( iparm[IPARM_AMALGAMATION_LVLCBLK] );
    error += iparm_reordering_split_check_value( iparm[IPARM_REORDERING_SPLIT] );
    error += iparm_reordering_stop_check_value( iparm[IPARM_REORDERING_STOP] );
    error += iparm_splitting_strategy_check_value( iparm[IPARM_SPLITTING_STRATEGY] );
    error += iparm_splitting_levels_projections_check_value( iparm[IPARM_SPLITTING_LEVELS_PROJECTIONS] );
    error += iparm_splitting_levels_kway_check_value( iparm[IPARM_SPLITTING_LEVELS_KWAY] );
    error += iparm_splitting_projections_depth_check_value( iparm[IPARM_SPLITTING_PROJECTIONS_DEPTH] );
    error += iparm_splitting_projections_distance_check_value( iparm[IPARM_SPLITTING_PROJECTIONS_DISTANCE] );
    error += iparm_splitting_projections_width_check_value( iparm[IPARM_SPLITTING_PROJECTIONS_WIDTH] );
    error += iparm_min_blocksize_check_value( iparm[IPARM_MIN_BLOCKSIZE] );
    error += iparm_max_blocksize_check_value( iparm[IPARM_MAX_BLOCKSIZE] );
    error += iparm_tasks2d_level_check_value( iparm[IPARM_TASKS2D_LEVEL] );
    error += iparm_tasks2d_width_check_value( iparm[IPARM_TASKS2D_WIDTH] );
    error += iparm_allcand_check_value( iparm[IPARM_ALLCAND] );
    error += iparm_incomplete_check_value( iparm[IPARM_INCOMPLETE] );
    error += iparm_level_of_fill_check_value( iparm[IPARM_LEVEL_OF_FILL] );
    error += iparm_factorization_check_value( iparm[IPARM_FACTORIZATION] );
    error += iparm_facto_look_side_check_value( iparm[IPARM_FACTO_LOOK_SIDE] );
    error += iparm_free_cscuser_check_value( iparm[IPARM_FREE_CSCUSER] );
    error += iparm_schur_fact_mode_check_value( iparm[IPARM_SCHUR_FACT_MODE] );
    error += iparm_transpose_solve_check_value( iparm[IPARM_TRANSPOSE_SOLVE] );
    error += iparm_schur_solv_mode_check_value( iparm[IPARM_SCHUR_SOLV_MODE] );
    error += iparm_applyperm_ws_check_value( iparm[IPARM_APPLYPERM_WS] );
    error += iparm_refinement_check_value( iparm[IPARM_REFINEMENT] );
    error += iparm_itermax_check_value( iparm[IPARM_ITERMAX] );
    error += iparm_gmres_im_check_value( iparm[IPARM_GMRES_IM] );
    error += iparm_scheduler_check_value( iparm[IPARM_SCHEDULER] );
    error += iparm_thread_nbr_check_value( iparm[IPARM_THREAD_NBR] );
    error += iparm_autosplit_comm_check_value( iparm[IPARM_AUTOSPLIT_COMM] );
    error += iparm_gpu_nbr_check_value( iparm[IPARM_GPU_NBR] );
    error += iparm_gpu_memory_percentage_check_value( iparm[IPARM_GPU_MEMORY_PERCENTAGE] );
    error += iparm_gpu_memory_block_size_check_value( iparm[IPARM_GPU_MEMORY_BLOCK_SIZE] );
    error += iparm_global_allocation_check_value( iparm[IPARM_GLOBAL_ALLOCATION] );
    error += iparm_compress_min_width_check_value( iparm[IPARM_COMPRESS_MIN_WIDTH] );
    error += iparm_compress_min_height_check_value( iparm[IPARM_COMPRESS_MIN_HEIGHT] );
    error += iparm_compress_when_check_value( iparm[IPARM_COMPRESS_WHEN] );
    error += iparm_compress_method_check_value( iparm[IPARM_COMPRESS_METHOD] );
    error += iparm_compress_ortho_check_value( iparm[IPARM_COMPRESS_ORTHO] );
    error += iparm_compress_reltol_check_value( iparm[IPARM_COMPRESS_RELTOL] );
    error += iparm_compress_preselect_check_value( iparm[IPARM_COMPRESS_PRESELECT] );
    error += iparm_compress_iluk_check_value( iparm[IPARM_COMPRESS_ILUK] );
    error += iparm_mixed_check_value( iparm[IPARM_MIXED] );
    error += iparm_ftz_check_value( iparm[IPARM_FTZ] );
    error += iparm_modify_parameter_check_value( iparm[IPARM_MODIFY_PARAMETER] );
    error += iparm_start_task_check_value( iparm[IPARM_START_TASK] );
    error += iparm_float_check_value( iparm[IPARM_FLOAT] );
    error += iparm_mtx_type_check_value( iparm[IPARM_MTX_TYPE] );
    error += iparm_dof_nbr_check_value( iparm[IPARM_DOF_NBR] );
    return error;
}

/**
 *******************************************************************************
 *
 * @brief Check the values of the dparm.
 *
 *******************************************************************************
 *
 * @param[in] dparm
 *          The dparm options array.
 *
 *******************************************************************************
 *
 * @return The amount of incorrect values in the dparm array.
 *
 *******************************************************************************/
int
dparm_check_values( const double *dparm )
{
    int error = 0;
    error += dparm_epsilon_refinement_check_value( dparm[DPARM_EPSILON_REFINEMENT] );
    error += dparm_epsilon_magn_ctrl_check_value( dparm[DPARM_EPSILON_MAGN_CTRL] );
    error += dparm_compress_tolerance_check_value( dparm[DPARM_COMPRESS_TOLERANCE] );
    error += dparm_compress_min_ratio_check_value( dparm[DPARM_COMPRESS_MIN_RATIO] );
    return error;
}
