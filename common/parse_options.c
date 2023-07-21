/**
 *
 * @file parse_options.c
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
#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */
#include <string.h>

/**
 *******************************************************************************
 *
 * @ingroup pastix_common
 *
 * @brief Parse iparm keywords to return its associated index in the iparm array
 *
 * This function converts the string only for input parameters, output
 * parameters are not handled.
 *
 *******************************************************************************
 *
 * @param[in] iparm
 *          The iparm string to convert to enum.
 *
 *******************************************************************************
 *
 * @retval The index of the iparm in the array.
 * @retval -1 if the string is not an iparm parameter.
 *
 *******************************************************************************/
pastix_iparm_t
parse_iparm( const char *iparm )
{
    if(0 == strcasecmp("iparm_verbose",                        iparm)) { return IPARM_VERBOSE; }
    if(0 == strcasecmp("iparm_io_strategy",                    iparm)) { return IPARM_IO_STRATEGY; }

    if(0 == strcasecmp("iparm_produce_stats",                  iparm)) { return IPARM_PRODUCE_STATS; }
    if(0 == strcasecmp("iparm_trace",                          iparm)) { return IPARM_TRACE; }

    if(0 == strcasecmp("iparm_mc64",                           iparm)) { return IPARM_MC64; }

    if(0 == strcasecmp("iparm_ordering",                       iparm)) { return IPARM_ORDERING; }
    if(0 == strcasecmp("iparm_ordering_default",               iparm)) { return IPARM_ORDERING_DEFAULT; }

    if(0 == strcasecmp("iparm_scotch_mt",                      iparm)) { return IPARM_SCOTCH_MT; }
    if(0 == strcasecmp("iparm_scotch_switch_level",            iparm)) { return IPARM_SCOTCH_SWITCH_LEVEL; }
    if(0 == strcasecmp("iparm_scotch_cmin",                    iparm)) { return IPARM_SCOTCH_CMIN; }
    if(0 == strcasecmp("iparm_scotch_cmax",                    iparm)) { return IPARM_SCOTCH_CMAX; }
    if(0 == strcasecmp("iparm_scotch_frat",                    iparm)) { return IPARM_SCOTCH_FRAT; }

    if(0 == strcasecmp("iparm_metis_ctype",                    iparm)) { return IPARM_METIS_CTYPE; }
    if(0 == strcasecmp("iparm_metis_rtype",                    iparm)) { return IPARM_METIS_RTYPE; }
    if(0 == strcasecmp("iparm_metis_no2hop",                   iparm)) { return IPARM_METIS_NO2HOP; }
    if(0 == strcasecmp("iparm_metis_nseps",                    iparm)) { return IPARM_METIS_NSEPS; }
    if(0 == strcasecmp("iparm_metis_niter",                    iparm)) { return IPARM_METIS_NITER; }
    if(0 == strcasecmp("iparm_metis_ufactor",                  iparm)) { return IPARM_METIS_UFACTOR; }
    if(0 == strcasecmp("iparm_metis_compress",                 iparm)) { return IPARM_METIS_COMPRESS; }
    if(0 == strcasecmp("iparm_metis_ccorder",                  iparm)) { return IPARM_METIS_CCORDER; }
    if(0 == strcasecmp("iparm_metis_pfactor",                  iparm)) { return IPARM_METIS_PFACTOR; }
    if(0 == strcasecmp("iparm_metis_seed",                     iparm)) { return IPARM_METIS_SEED; }
    if(0 == strcasecmp("iparm_metis_dbglvl",                   iparm)) { return IPARM_METIS_DBGLVL; }

    if(0 == strcasecmp("iparm_amalgamation_lvlblas",           iparm)) { return IPARM_AMALGAMATION_LVLBLAS; }
    if(0 == strcasecmp("iparm_amalgamation_lvlcblk",           iparm)) { return IPARM_AMALGAMATION_LVLCBLK; }

    if(0 == strcasecmp("iparm_reordering_split",               iparm)) { return IPARM_REORDERING_SPLIT; }
    if(0 == strcasecmp("iparm_reordering_stop",                iparm)) { return IPARM_REORDERING_STOP; }
    if(0 == strcasecmp("iparm_splitting_strategy",             iparm)) { return IPARM_SPLITTING_STRATEGY; }
    if(0 == strcasecmp("iparm_splitting_levels_projections",   iparm)) { return IPARM_SPLITTING_LEVELS_PROJECTIONS; }
    if(0 == strcasecmp("iparm_splitting_levels_kway",          iparm)) { return IPARM_SPLITTING_LEVELS_KWAY; }
    if(0 == strcasecmp("iparm_splitting_projections_depth",    iparm)) { return IPARM_SPLITTING_PROJECTIONS_DEPTH; }
    if(0 == strcasecmp("iparm_splitting_projections_distance", iparm)) { return IPARM_SPLITTING_PROJECTIONS_DISTANCE; }
    if(0 == strcasecmp("iparm_splitting_projections_width",    iparm)) { return IPARM_SPLITTING_PROJECTIONS_WIDTH; }

    if(0 == strcasecmp("iparm_min_blocksize",                  iparm)) { return IPARM_MIN_BLOCKSIZE; }
    if(0 == strcasecmp("iparm_max_blocksize",                  iparm)) { return IPARM_MAX_BLOCKSIZE; }
    if(0 == strcasecmp("iparm_tasks2d_level",                  iparm)) { return IPARM_TASKS2D_LEVEL; }
    if(0 == strcasecmp("iparm_tasks2d_width",                  iparm)) { return IPARM_TASKS2D_WIDTH; }
    if(0 == strcasecmp("iparm_allcand",                        iparm)) { return IPARM_ALLCAND; }

    if(0 == strcasecmp("iparm_incomplete",                     iparm)) { return IPARM_INCOMPLETE; }
    if(0 == strcasecmp("iparm_level_of_fill",                  iparm)) { return IPARM_LEVEL_OF_FILL; }

    if(0 == strcasecmp("iparm_factorization",                  iparm)) { return IPARM_FACTORIZATION; }
    if(0 == strcasecmp("iparm_facto_look_side",                iparm)) { return IPARM_FACTO_LOOK_SIDE; }
    if(0 == strcasecmp("iparm_free_cscuser",                   iparm)) { return IPARM_FREE_CSCUSER; }
    if(0 == strcasecmp("iparm_schur_fact_mode",                iparm)) { return IPARM_SCHUR_FACT_MODE; }

    if(0 == strcasecmp("iparm_transpose_solve",                iparm)) { return IPARM_TRANSPOSE_SOLVE; }
    if(0 == strcasecmp("iparm_schur_solv_mode",                iparm)) { return IPARM_SCHUR_SOLV_MODE; }
    if(0 == strcasecmp("iparm_applyperm_ws",                   iparm)) { return IPARM_APPLYPERM_WS; }

    if(0 == strcasecmp("iparm_refinement",                     iparm)) { return IPARM_REFINEMENT; }
    if(0 == strcasecmp("iparm_itermax",                        iparm)) { return IPARM_ITERMAX; }
    if(0 == strcasecmp("iparm_gmres_im",                       iparm)) { return IPARM_GMRES_IM; }

    if(0 == strcasecmp("iparm_scheduler",                      iparm)) { return IPARM_SCHEDULER; }
    if(0 == strcasecmp("iparm_thread_nbr",                     iparm)) { return IPARM_THREAD_NBR; }
    if(0 == strcasecmp("iparm_autosplit_comm",                 iparm)) { return IPARM_AUTOSPLIT_COMM; }

    if(0 == strcasecmp("iparm_gpu_nbr",                        iparm)) { return IPARM_GPU_NBR; }
    if(0 == strcasecmp("iparm_gpu_memory_percentage",          iparm)) { return IPARM_GPU_MEMORY_PERCENTAGE; }
    if(0 == strcasecmp("iparm_gpu_memory_block_size",          iparm)) { return IPARM_GPU_MEMORY_BLOCK_SIZE; }
    if(0 == strcasecmp("iparm_global_allocation",              iparm)) { return IPARM_GLOBAL_ALLOCATION; }

    if(0 == strcasecmp("iparm_compress_min_width",             iparm)) { return IPARM_COMPRESS_MIN_WIDTH; }
    if(0 == strcasecmp("iparm_compress_min_height",            iparm)) { return IPARM_COMPRESS_MIN_HEIGHT; }
    if(0 == strcasecmp("iparm_compress_when",                  iparm)) { return IPARM_COMPRESS_WHEN; }
    if(0 == strcasecmp("iparm_compress_method",                iparm)) { return IPARM_COMPRESS_METHOD; }
    if(0 == strcasecmp("iparm_compress_ortho",                 iparm)) { return IPARM_COMPRESS_ORTHO; }
    if(0 == strcasecmp("iparm_compress_reltol",                iparm)) { return IPARM_COMPRESS_RELTOL; }
    if(0 == strcasecmp("iparm_compress_preselect",             iparm)) { return IPARM_COMPRESS_PRESELECT; }
    if(0 == strcasecmp("iparm_compress_iluk",                  iparm)) { return IPARM_COMPRESS_ILUK; }

    if(0 == strcasecmp("iparm_mixed",                          iparm)) { return IPARM_MIXED; }
    if(0 == strcasecmp("iparm_ftz",                            iparm)) { return IPARM_FTZ; }

    if(0 == strcasecmp("iparm_modify_parameter",               iparm)) { return IPARM_MODIFY_PARAMETER; }
    if(0 == strcasecmp("iparm_start_task",                     iparm)) { return IPARM_START_TASK; }
    if(0 == strcasecmp("iparm_end_task",                       iparm)) { return IPARM_END_TASK; }
    if(0 == strcasecmp("iparm_float",                          iparm)) { return IPARM_FLOAT; }
    if(0 == strcasecmp("iparm_mtx_type",                       iparm)) { return IPARM_MTX_TYPE; }
    if(0 == strcasecmp("iparm_dof_nbr",                        iparm)) { return IPARM_DOF_NBR; }

    return -1;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_common
 *
 * @brief Parse dparm keywords to return its associated index in the dparm array
 *
 * This function converts the string only for input parameters, output
 * parameters are not handled.
 *
 *******************************************************************************
 *
 * @param[in] dparm
 *          The dparm string to convert to enum.
 *
 *******************************************************************************
 *
 * @retval The index of the dparm in the array.
 * @retval -1 if the string is not a dparm parameter.
 *
 *******************************************************************************/
pastix_dparm_t
parse_dparm( const char *dparm )
{
    if(0 == strcasecmp("dparm_epsilon_refinement", dparm)) { return DPARM_EPSILON_REFINEMENT; }
    if(0 == strcasecmp("dparm_epsilon_magn_ctrl",  dparm)) { return DPARM_EPSILON_MAGN_CTRL; }
    if(0 == strcasecmp("dparm_compress_tolerance", dparm)) { return DPARM_COMPRESS_TOLERANCE; }
    if(0 == strcasecmp("dparm_compress_min_ratio", dparm)) { return DPARM_COMPRESS_MIN_RATIO; }

    return -1;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_common
 *
 * @brief Parse enum values for iparm array, and return the enum value
 * associated to it.
 *
 *******************************************************************************
 *
 * @param[in] string
 *          The enum name to convert to its value
 *
 *******************************************************************************
 *
 * @retval The value if the enum associated to the string
 * @retval -1 if the string is not an enum in the pastix API.
 *
 *******************************************************************************/
int
parse_enums( const char *string )
{
    if(0 == strcasecmp("pastixverbosenot", string)) { return PastixVerboseNot; }
    if(0 == strcasecmp("pastixverboseno",  string)) { return PastixVerboseNo; }
    if(0 == strcasecmp("pastixverboseyes", string)) { return PastixVerboseYes; }

    if(0 == strcasecmp("pastixiono",        string)) { return PastixIONo; }
    if(0 == strcasecmp("pastixioload",      string)) { return PastixIOLoad; }
    if(0 == strcasecmp("pastixiosave",      string)) { return PastixIOSave; }
    if(0 == strcasecmp("pastixioloadgraph", string)) { return PastixIOLoadGraph; }
    if(0 == strcasecmp("pastixiosavegraph", string)) { return PastixIOSaveGraph; }
    if(0 == strcasecmp("pastixioloadcsc",   string)) { return PastixIOLoadCSC; }
    if(0 == strcasecmp("pastixiosavecsc",   string)) { return PastixIOSaveCSC; }

    if(0 == strcasecmp("pastixtracenumfact", string)) { return PastixTraceNumfact; }
    if(0 == strcasecmp("pastixtracesolve",   string)) { return PastixTraceSolve; }

    if(0 == strcasecmp("pastixorderscotch",   string)) { return PastixOrderScotch; }
    if(0 == strcasecmp("pastixordermetis",    string)) { return PastixOrderMetis; }
    if(0 == strcasecmp("pastixorderpersonal", string)) { return PastixOrderPersonal; }
    if(0 == strcasecmp("pastixorderptscotch", string)) { return PastixOrderPtScotch; }
    if(0 == strcasecmp("pastixorderparmetis", string)) { return PastixOrderParMetis; }

    if(0 == strcasecmp("pastixsplitnot",             string)) { return PastixSplitNot; }
    if(0 == strcasecmp("pastixsplitkway",            string)) { return PastixSplitKway; }
    if(0 == strcasecmp("pastixsplitkwayprojections", string)) { return PastixSplitKwayProjections; }

    if(0 == strcasecmp("pastixfactpotrf", string)) { return PastixFactPOTRF; }
    if(0 == strcasecmp("pastixfactsytrf", string)) { return PastixFactSYTRF; }
    if(0 == strcasecmp("pastixfactgetrf", string)) { return PastixFactGETRF; }
    if(0 == strcasecmp("pastixfactpxtrf", string)) { return PastixFactPXTRF; }
    if(0 == strcasecmp("pastixfacthetrf", string)) { return PastixFactHETRF; }
    if(0 == strcasecmp("pastixfactllh",   string)) { return PastixFactLLH; }
    if(0 == strcasecmp("pastixfactldlt",  string)) { return PastixFactLDLT; }
    if(0 == strcasecmp("pastixfactlu",    string)) { return PastixFactLU; }
    if(0 == strcasecmp("pastixfactllt",   string)) { return PastixFactLLT; }
    if(0 == strcasecmp("pastixfactldlh",  string)) { return PastixFactLDLH; }

    if(0 == strcasecmp("pastixfactmodelocal", string)) { return PastixFactModeLocal; }
    if(0 == strcasecmp("pastixfactmodeschur", string)) { return PastixFactModeSchur; }
    if(0 == strcasecmp("pastixfactmodeboth",  string)) { return PastixFactModeBoth; }

    if(0 == strcasecmp("pastixnotrans",   string)) { return PastixNoTrans; }
    if(0 == strcasecmp("pastixtrans",     string)) { return PastixTrans; }
    if(0 == strcasecmp("pastixconjtrans", string)) { return PastixConjTrans; }

    if(0 == strcasecmp("pastixsolvmodelocal",     string)) { return PastixSolvModeLocal; }
    if(0 == strcasecmp("pastixsolvmodeinterface", string)) { return PastixSolvModeInterface; }
    if(0 == strcasecmp("pastixsolvmodeschur",     string)) { return PastixSolvModeSchur; }

    if(0 == strcasecmp("pastixrefinegmres",    string)) { return PastixRefineGMRES; }
    if(0 == strcasecmp("pastixrefinecg",       string)) { return PastixRefineCG; }
    if(0 == strcasecmp("pastixrefinesr",       string)) { return PastixRefineSR; }
    if(0 == strcasecmp("pastixrefinebicgstab", string)) { return PastixRefineBiCGSTAB; }

    if(0 == strcasecmp("pastixschedsequential", string)) { return PastixSchedSequential; }
    if(0 == strcasecmp("pastixschedstatic",     string)) { return PastixSchedStatic; }
    if(0 == strcasecmp("pastixschedparsec",     string)) { return PastixSchedParsec; }
    if(0 == strcasecmp("pastixschedstarpu",     string)) { return PastixSchedStarPU; }
    if(0 == strcasecmp("pastixscheddynamic",    string)) { return PastixSchedDynamic; }

    if(0 == strcasecmp("pastixcompressnever",      string)) { return PastixCompressNever; }
    if(0 == strcasecmp("pastixcompresswhenbegin",  string)) { return PastixCompressWhenBegin; }
    if(0 == strcasecmp("pastixcompresswhenend",    string)) { return PastixCompressWhenEnd; }
    if(0 == strcasecmp("pastixcompresswhenduring", string)) { return PastixCompressWhenDuring; }

    if(0 == strcasecmp("pastixcompressmethodsvd",   string)) { return PastixCompressMethodSVD; }
    if(0 == strcasecmp("pastixcompressmethodpqrcp", string)) { return PastixCompressMethodPQRCP; }
    if(0 == strcasecmp("pastixcompressmethodrqrcp", string)) { return PastixCompressMethodRQRCP; }
    if(0 == strcasecmp("pastixcompressmethodtqrcp", string)) { return PastixCompressMethodTQRCP; }
    if(0 == strcasecmp("pastixcompressmethodrqrrt", string)) { return PastixCompressMethodRQRRT; }
    if(0 == strcasecmp("pastixcompressmethodnbr",   string)) { return PastixCompressMethodNbr; }

    if(0 == strcasecmp("pastixcompressorthocgs",       string)) { return PastixCompressOrthoCGS; }
    if(0 == strcasecmp("pastixcompressorthoqr",        string)) { return PastixCompressOrthoQR; }
    if(0 == strcasecmp("pastixcompressorthopartialqr", string)) { return PastixCompressOrthoPartialQR; }

    if(0 == strcasecmp("pastixtaskinit",     string)) { return PastixTaskInit; }
    if(0 == strcasecmp("pastixtaskordering", string)) { return PastixTaskOrdering; }
    if(0 == strcasecmp("pastixtasksymbfact", string)) { return PastixTaskSymbfact; }
    if(0 == strcasecmp("pastixtaskanalyze",  string)) { return PastixTaskAnalyze; }
    if(0 == strcasecmp("pastixtasknumfact",  string)) { return PastixTaskNumfact; }
    if(0 == strcasecmp("pastixtasksolve",    string)) { return PastixTaskSolve; }
    if(0 == strcasecmp("pastixtaskrefine",   string)) { return PastixTaskRefine; }
    if(0 == strcasecmp("pastixtaskclean",    string)) { return PastixTaskClean; }

    if(0 == strcasecmp("pastixpattern",   string)) { return PastixPattern; }
    if(0 == strcasecmp("pastixfloat",     string)) { return PastixFloat; }
    if(0 == strcasecmp("pastixdouble",    string)) { return PastixDouble; }
    if(0 == strcasecmp("pastixcomplex32", string)) { return PastixComplex32; }
    if(0 == strcasecmp("pastixcomplex64", string)) { return PastixComplex64; }

    /* If the value is directly given without string */
    {
        int value;
        if ( sscanf( string, "%d", &value ) != 1 ) {
            return -1;
        }
        else {
            return value;
        }
    }
}

const char*
pastix_task_getstr( pastix_task_t value )
{
    switch( value ) {
    case PastixTaskInit:
        return "PastixTaskInit";
    case PastixTaskOrdering:
        return "PastixTaskOrdering";
    case PastixTaskSymbfact:
        return "PastixTaskSymbfact";
    case PastixTaskAnalyze:
        return "PastixTaskAnalyze";
    case PastixTaskNumfact:
        return "PastixTaskNumfact";
    case PastixTaskSolve:
        return "PastixTaskSolve";
    case PastixTaskRefine:
        return "PastixTaskRefine";
    case PastixTaskClean:
        return "PastixTaskClean";
    default :
        return "Bad task given";
    }
}

const char*
pastix_trace_getstr( pastix_trace_t value )
{
    switch( value ) {
    case PastixTraceNumfact:
        return "PastixTraceNumfact";
    case PastixTraceSolve:
        return "PastixTraceSolve";
    default :
        return "Bad trace given";
    }
}

const char*
pastix_verbose_getstr( pastix_verbose_t value )
{
    switch( value ) {
    case PastixVerboseNot:
        return "PastixVerboseNot";
    case PastixVerboseNo:
        return "PastixVerboseNo";
    case PastixVerboseYes:
        return "PastixVerboseYes";
    default :
        return "Bad verbose given";
    }
}

const char*
pastix_io_getstr( pastix_io_t value )
{
    switch( value ) {
    case PastixIONo:
        return "PastixIONo";
    case PastixIOLoad:
        return "PastixIOLoad";
    case PastixIOSave:
        return "PastixIOSave";
    case PastixIOLoadGraph:
        return "PastixIOLoadGraph";
    case PastixIOSaveGraph:
        return "PastixIOSaveGraph";
    case PastixIOLoadCSC:
        return "PastixIOLoadCSC";
    case PastixIOSaveCSC:
        return "PastixIOSaveCSC";
    default :
        return "Bad io given";
    }
}

const char*
pastix_fact_mode_getstr( pastix_fact_mode_t value )
{
    switch( value ) {
    case PastixFactModeLocal:
        return "PastixFactModeLocal";
    case PastixFactModeSchur:
        return "PastixFactModeSchur";
    case PastixFactModeBoth:
        return "PastixFactModeBoth";
    default :
        return "Bad fact_mode given";
    }
}

const char*
pastix_solv_mode_getstr( pastix_solv_mode_t value )
{
    switch( value ) {
    case PastixSolvModeLocal:
        return "PastixSolvModeLocal";
    case PastixSolvModeInterface:
        return "PastixSolvModeInterface";
    case PastixSolvModeSchur:
        return "PastixSolvModeSchur";
    default :
        return "Bad solv_mode given";
    }
}

const char*
pastix_refine_getstr( pastix_refine_t value )
{
    switch( value ) {
    case PastixRefineGMRES:
        return "PastixRefineGMRES";
    case PastixRefineCG:
        return "PastixRefineCG";
    case PastixRefineSR:
        return "PastixRefineSR";
    case PastixRefineBiCGSTAB:
        return "PastixRefineBiCGSTAB";
    default :
        return "Bad refine given";
    }
}

const char*
pastix_coeftype_getstr( pastix_coeftype_t value )
{
    switch( value ) {
    case PastixPattern:
        return "PastixPattern";
    case PastixFloat:
        return "PastixFloat";
    case PastixDouble:
        return "PastixDouble";
    case PastixComplex32:
        return "PastixComplex32";
    case PastixComplex64:
        return "PastixComplex64";
    default :
        return "Bad coeftype given";
    }
}

const char*
pastix_factotype_getstr( pastix_factotype_t value )
{
    switch( value ) {
    case PastixFactLLH:
        return "PastixFactLLH";
    case PastixFactLDLT:
        return "PastixFactLDLT";
    case PastixFactLU:
        return "PastixFactLU";
    case PastixFactLLT:
        return "PastixFactLLT";
    case PastixFactLDLH:
        return "PastixFactLDLH";
    default :
        return "Bad factotype given";
    }
}

const char*
pastix_factolookside_getstr( pastix_factolookside_t value )
{
    switch( value ) {
    case PastixFactLeftLooking:
        return "PastixFactLeftLooking";
    case PastixFactRightLooking:
        return "PastixFactRightLooking";
    default :
        return "Bad facto looking side given";
    }
}

const char*
pastix_scheduler_getstr( pastix_scheduler_t value )
{
    switch( value ) {
    case PastixSchedSequential:
        return "PastixSchedSequential";
    case PastixSchedStatic:
        return "PastixSchedStatic";
    case PastixSchedParsec:
        return "PastixSchedParsec";
    case PastixSchedStarPU:
        return "PastixSchedStarPU";
    case PastixSchedDynamic:
        return "PastixSchedDynamic";
    default :
        return "Bad scheduler given";
    }
}

const char*
pastix_ordering_getstr( pastix_ordering_t value )
{
    switch( value ) {
    case PastixOrderScotch:
        return "PastixOrderScotch";
    case PastixOrderMetis:
        return "PastixOrderMetis";
    case PastixOrderPersonal:
        return "PastixOrderPersonal";
    case PastixOrderPtScotch:
        return "PastixOrderPtScotch";
    case PastixOrderParMetis:
        return "PastixOrderParMetis";
    default :
        return "Bad ordering given";
    }
}

const char*
pastix_mpithreadmode_getstr( pastix_mpithreadmode_t value )
{
    switch( value ) {
    case PastixMpiNone:
        return "PastixMpiNone";
    case PastixMpiThreadSingle:
        return "PastixMpiThreadSingle";
    case PastixMpiThreadFunneled:
        return "PastixMpiThreadFunneled";
    case PastixMpiThreadSerialized:
        return "PastixMpiThreadSerialized";
    case PastixMpiThreadMultiple:
        return "PastixMpiThreadMultiple";
    default :
        return "Bad mpithreadmode given";
    }
}

const char*
pastix_error_getstr( pastix_error_t value )
{
    switch( value ) {
    case PASTIX_SUCCESS:
        return "PASTIX_SUCCESS";
    case PASTIX_ERR_UNKNOWN:
        return "PASTIX_ERR_UNKNOWN";
    case PASTIX_ERR_ALLOC:
        return "PASTIX_ERR_ALLOC";
    case PASTIX_ERR_NOTIMPLEMENTED:
        return "PASTIX_ERR_NOTIMPLEMENTED";
    case PASTIX_ERR_OUTOFMEMORY:
        return "PASTIX_ERR_OUTOFMEMORY";
    case PASTIX_ERR_THREAD:
        return "PASTIX_ERR_THREAD";
    case PASTIX_ERR_INTERNAL:
        return "PASTIX_ERR_INTERNAL";
    case PASTIX_ERR_BADPARAMETER:
        return "PASTIX_ERR_BADPARAMETER";
    case PASTIX_ERR_FILE:
        return "PASTIX_ERR_FILE";
    case PASTIX_ERR_INTEGER_TYPE:
        return "PASTIX_ERR_INTEGER_TYPE";
    case PASTIX_ERR_IO:
        return "PASTIX_ERR_IO";
    case PASTIX_ERR_MPI:
        return "PASTIX_ERR_MPI";
    default :
        return "Bad error given";
    }
}

const char*
pastix_compress_when_getstr( pastix_compress_when_t value )
{
    switch( value ) {
    case PastixCompressNever:
        return "PastixCompressNever";
    case PastixCompressWhenBegin:
        return "PastixCompressWhenBegin";
    case PastixCompressWhenEnd:
        return "PastixCompressWhenEnd";
    case PastixCompressWhenDuring:
        return "PastixCompressWhenDuring";
    default :
        return "Bad compress_when given";
    }
}

const char*
pastix_compress_method_getstr( pastix_compress_method_t value )
{
    switch( value ) {
    case PastixCompressMethodSVD:
        return "PastixCompressMethodSVD";
    case PastixCompressMethodPQRCP:
        return "PastixCompressMethodPQRCP";
    case PastixCompressMethodRQRCP:
        return "PastixCompressMethodRQRCP";
    case PastixCompressMethodTQRCP:
        return "PastixCompressMethodTQRCP";
    case PastixCompressMethodRQRRT:
        return "PastixCompressMethodRQRRT";
    case PastixCompressMethodNbr:
        return "PastixCompressMethodNbr";
    default :
        return "Bad compress_method given";
    }
}

const char*
pastix_compress_ortho_getstr( pastix_compress_ortho_t value )
{
    switch( value ) {
    case PastixCompressOrthoCGS:
        return "PastixCompressOrthoCGS";
    case PastixCompressOrthoQR:
        return "PastixCompressOrthoQR";
    case PastixCompressOrthoPartialQR:
        return "PastixCompressOrthoPartialQR";
    default :
        return "Bad compress_ortho given";
    }
}

const char*
pastix_split_getstr( pastix_split_t value )
{
    switch( value ) {
    case PastixSplitNot:
        return "PastixSplitNot";
    case PastixSplitKway:
        return "PastixSplitKway";
    case PastixSplitKwayProjections:
        return "PastixSplitKwayProjections";
    default :
        return "Bad split given";
    }
}

const char*
pastix_layout_getstr( pastix_layout_t value )
{
    switch( value ) {
    case PastixRowMajor:
        return "PastixRowMajor";
    case PastixColMajor:
        return "PastixColMajor";
    default :
        return "Bad layout given";
    }
}

const char*
pastix_trans_getstr( pastix_trans_t value )
{
    switch( value ) {
    case PastixNoTrans:
        return "PastixNoTrans";
    case PastixTrans:
        return "PastixTrans";
    case PastixConjTrans:
        return "PastixConjTrans";
    default :
        return "Bad trans given";
    }
}

const char*
pastix_mtxtype_getstr( pastix_mtxtype_t value )
{
    switch( value ) {
    case PastixGeneral:
        return "PastixGeneral";
    case PastixSymmetric:
        return "PastixSymmetric";
    case PastixHermitian:
        return "PastixHermitian";
    default :
        return "Bad mtxtype given";
    }
}

const char*
pastix_uplo_getstr( pastix_uplo_t value )
{
    switch( value ) {
    case PastixUpper:
        return "PastixUpper";
    case PastixLower:
        return "PastixLower";
    case PastixUpperLower:
        return "PastixUpperLower";
    default :
        return "Bad uplo given";
    }
}

const char*
pastix_coefside_getstr( pastix_coefside_t value )
{
    switch( value ) {
    case PastixLCoef:
        return "PastixLCoef";
    case PastixUCoef:
        return "PastixUCoef";
    case PastixLUCoef:
        return "PastixLUCoef";
    default :
        return "Bad coefside given";
    }
}

const char*
pastix_diag_getstr( pastix_diag_t value )
{
    switch( value ) {
    case PastixNonUnit:
        return "PastixNonUnit";
    case PastixUnit:
        return "PastixUnit";
    default :
        return "Bad diag given";
    }
}

const char*
pastix_side_getstr( pastix_side_t value )
{
    switch( value ) {
    case PastixLeft:
        return "PastixLeft";
    case PastixRight:
        return "PastixRight";
    default :
        return "Bad side given";
    }
}

const char*
pastix_normtype_getstr( pastix_normtype_t value )
{
    switch( value ) {
    case PastixOneNorm:
        return "PastixOneNorm";
    case PastixFrobeniusNorm:
        return "PastixFrobeniusNorm";
    case PastixInfNorm:
        return "PastixInfNorm";
    case PastixMaxNorm:
        return "PastixMaxNorm";
    default :
        return "Bad normtype given";
    }
}

const char*
pastix_dir_getstr( pastix_dir_t value )
{
    switch( value ) {
    case PastixDirForward:
        return "PastixDirForward";
    case PastixDirBackward:
        return "PastixDirBackward";
    default :
        return "Bad dir given";
    }
}

/**
 *******************************************************************************
 *
 * @brief Dump the iparm an dparm parameters in the CSV file.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The main data structure.
 *
 * @param[inout] csv
 *          The csv file that will contain the dumped datas.
 *
 *******************************************************************************/
void
pastix_param2csv( const pastix_data_t *pastix_data,
                        FILE          *csv )
{
    pastix_int_t *iparm = pastix_data->iparm;
    double       *dparm = pastix_data->dparm;

    fprintf( csv, "%s,%s\n",  "iparm_verbose",      pastix_verbose_getstr(iparm[IPARM_VERBOSE]) );
    fprintf( csv, "%s,%s\n",  "iparm_io_strategy",  pastix_io_getstr(iparm[IPARM_IO_STRATEGY]) );

    fprintf( csv, "%s,%ld\n", "iparm_nnzeros",             (long)iparm[IPARM_NNZEROS] );
    fprintf( csv, "%s,%ld\n", "iparm_nnzeros_block_local", (long)iparm[IPARM_NNZEROS_BLOCK_LOCAL] );
    fprintf( csv, "%s,%ld\n", "iparm_allocated_terms",     (long)iparm[IPARM_ALLOCATED_TERMS] );
    fprintf( csv, "%s,%ld\n", "iparm_produce_stats",       (long)iparm[IPARM_PRODUCE_STATS] );
    fprintf( csv, "%s,%s\n",  "iparm_trace",                pastix_trace_getstr(iparm[IPARM_TRACE]) );

    fprintf( csv, "%s,%ld\n", "iparm_mc64", (long)iparm[IPARM_MC64] );

    fprintf( csv, "%s,%s\n",  "iparm_ordering",          pastix_ordering_getstr(iparm[IPARM_ORDERING]) );
    fprintf( csv, "%s,%ld\n", "iparm_ordering_default", (long)iparm[IPARM_ORDERING_DEFAULT] );

    fprintf( csv, "%s,%ld\n", "iparm_scotch_mt",           (long)iparm[IPARM_SCOTCH_MT] );
    fprintf( csv, "%s,%ld\n", "iparm_scotch_switch_level", (long)iparm[IPARM_SCOTCH_SWITCH_LEVEL] );
    fprintf( csv, "%s,%ld\n", "iparm_scotch_cmin",         (long)iparm[IPARM_SCOTCH_CMIN] );
    fprintf( csv, "%s,%ld\n", "iparm_scotch_cmax",         (long)iparm[IPARM_SCOTCH_CMAX] );
    fprintf( csv, "%s,%ld\n", "iparm_scotch_frat",         (long)iparm[IPARM_SCOTCH_FRAT] );

    fprintf( csv, "%s,%ld\n", "iparm_metis_ctype",    (long)iparm[IPARM_METIS_CTYPE] );
    fprintf( csv, "%s,%ld\n", "iparm_metis_rtype",    (long)iparm[IPARM_METIS_RTYPE] );
    fprintf( csv, "%s,%ld\n", "iparm_metis_no2hop",   (long)iparm[IPARM_METIS_NO2HOP] );
    fprintf( csv, "%s,%ld\n", "iparm_metis_nseps",    (long)iparm[IPARM_METIS_NSEPS] );
    fprintf( csv, "%s,%ld\n", "iparm_metis_niter",    (long)iparm[IPARM_METIS_NITER] );
    fprintf( csv, "%s,%ld\n", "iparm_metis_ufactor",  (long)iparm[IPARM_METIS_UFACTOR] );
    fprintf( csv, "%s,%ld\n", "iparm_metis_compress", (long)iparm[IPARM_METIS_COMPRESS] );
    fprintf( csv, "%s,%ld\n", "iparm_metis_ccorder",  (long)iparm[IPARM_METIS_CCORDER] );
    fprintf( csv, "%s,%ld\n", "iparm_metis_pfactor",  (long)iparm[IPARM_METIS_PFACTOR] );
    fprintf( csv, "%s,%ld\n", "iparm_metis_seed",     (long)iparm[IPARM_METIS_SEED] );
    fprintf( csv, "%s,%ld\n", "iparm_metis_dbglvl",   (long)iparm[IPARM_METIS_DBGLVL] );

    fprintf( csv, "%s,%ld\n", "iparm_amalgamation_lvlblas", (long)iparm[IPARM_AMALGAMATION_LVLBLAS] );
    fprintf( csv, "%s,%ld\n", "iparm_amalgamation_lvlcblk", (long)iparm[IPARM_AMALGAMATION_LVLCBLK] );

    fprintf( csv, "%s,%ld\n", "iparm_reordering_split",               (long)iparm[IPARM_REORDERING_SPLIT] );
    fprintf( csv, "%s,%ld\n", "iparm_reordering_stop",                (long)iparm[IPARM_REORDERING_STOP] );
    fprintf( csv, "%s,%s\n",  "iparm_splitting_strategy",              pastix_split_getstr(iparm[IPARM_SPLITTING_STRATEGY]) );
    fprintf( csv, "%s,%ld\n", "iparm_splitting_levels_projections",   (long)iparm[IPARM_SPLITTING_LEVELS_PROJECTIONS] );
    fprintf( csv, "%s,%ld\n", "iparm_splitting_levels_kway",          (long)iparm[IPARM_SPLITTING_LEVELS_KWAY] );
    fprintf( csv, "%s,%ld\n", "iparm_splitting_projections_depth",    (long)iparm[IPARM_SPLITTING_PROJECTIONS_DEPTH] );
    fprintf( csv, "%s,%ld\n", "iparm_splitting_projections_distance", (long)iparm[IPARM_SPLITTING_PROJECTIONS_DISTANCE] );
    fprintf( csv, "%s,%ld\n", "iparm_splitting_projections_width",    (long)iparm[IPARM_SPLITTING_PROJECTIONS_WIDTH] );

    fprintf( csv, "%s,%ld\n", "iparm_min_blocksize", (long)iparm[IPARM_MIN_BLOCKSIZE] );
    fprintf( csv, "%s,%ld\n", "iparm_max_blocksize", (long)iparm[IPARM_MAX_BLOCKSIZE] );
    fprintf( csv, "%s,%ld\n", "iparm_tasks2d_level", (long)iparm[IPARM_TASKS2D_LEVEL] );
    fprintf( csv, "%s,%ld\n", "iparm_tasks2d_width", (long)iparm[IPARM_TASKS2D_WIDTH] );
    fprintf( csv, "%s,%ld\n", "iparm_allcand",       (long)iparm[IPARM_ALLCAND] );

    fprintf( csv, "%s,%ld\n", "iparm_incomplete",    (long)iparm[IPARM_INCOMPLETE] );
    fprintf( csv, "%s,%ld\n", "iparm_level_of_fill", (long)iparm[IPARM_LEVEL_OF_FILL] );

    fprintf( csv, "%s,%s\n",  "iparm_factorization",    pastix_factotype_getstr(iparm[IPARM_FACTORIZATION]) );
    fprintf( csv, "%s,%s\n",  "iparm_facto_look_side",  pastix_factolookside_getstr(iparm[IPARM_FACTO_LOOK_SIDE]) );
    fprintf( csv, "%s,%ld\n", "iparm_static_pivoting", (long)iparm[IPARM_STATIC_PIVOTING] );
    fprintf( csv, "%s,%ld\n", "iparm_free_cscuser",    (long)iparm[IPARM_FREE_CSCUSER] );
    fprintf( csv, "%s,%s\n",  "iparm_schur_fact_mode",  pastix_fact_mode_getstr(iparm[IPARM_SCHUR_FACT_MODE]) );

    fprintf( csv, "%s,%s\n",  "iparm_transpose_solve",  pastix_trans_getstr(iparm[IPARM_TRANSPOSE_SOLVE]) );
    fprintf( csv, "%s,%s\n",  "iparm_schur_solv_mode",  pastix_solv_mode_getstr(iparm[IPARM_SCHUR_SOLV_MODE]) );
    fprintf( csv, "%s,%ld\n", "iparm_applyperm_ws",    (long)iparm[IPARM_APPLYPERM_WS] );

    fprintf( csv, "%s,%s\n",  "iparm_refinement",  pastix_refine_getstr(iparm[IPARM_REFINEMENT]) );
    fprintf( csv, "%s,%ld\n", "iparm_nbiter",     (long)iparm[IPARM_NBITER] );
    fprintf( csv, "%s,%ld\n", "iparm_itermax",    (long)iparm[IPARM_ITERMAX] );
    fprintf( csv, "%s,%ld\n", "iparm_gmres_im",   (long)iparm[IPARM_GMRES_IM] );

    fprintf( csv, "%s,%s\n",  "iparm_scheduler",       pastix_scheduler_getstr(iparm[IPARM_SCHEDULER]) );
    fprintf( csv, "%s,%ld\n", "iparm_thread_nbr",     (long)iparm[IPARM_THREAD_NBR] );
    fprintf( csv, "%s,%ld\n", "iparm_autosplit_comm", (long)iparm[IPARM_AUTOSPLIT_COMM] );

    fprintf( csv, "%s,%ld\n", "iparm_gpu_nbr",               (long)iparm[IPARM_GPU_NBR] );
    fprintf( csv, "%s,%ld\n", "iparm_gpu_memory_percentage", (long)iparm[IPARM_GPU_MEMORY_PERCENTAGE] );
    fprintf( csv, "%s,%ld\n", "iparm_gpu_memory_block_size", (long)iparm[IPARM_GPU_MEMORY_BLOCK_SIZE] );
    fprintf( csv, "%s,%ld\n", "iparm_global_allocation",     (long)iparm[IPARM_GLOBAL_ALLOCATION] );

    fprintf( csv, "%s,%ld\n", "iparm_compress_min_width",  (long)iparm[IPARM_COMPRESS_MIN_WIDTH] );
    fprintf( csv, "%s,%ld\n", "iparm_compress_min_height", (long)iparm[IPARM_COMPRESS_MIN_HEIGHT] );
    fprintf( csv, "%s,%s\n",  "iparm_compress_when",        pastix_compress_when_getstr(iparm[IPARM_COMPRESS_WHEN]) );
    fprintf( csv, "%s,%s\n",  "iparm_compress_method",      pastix_compress_method_getstr(iparm[IPARM_COMPRESS_METHOD]) );
    fprintf( csv, "%s,%s\n",  "iparm_compress_ortho",       pastix_compress_ortho_getstr(iparm[IPARM_COMPRESS_ORTHO]) );
    fprintf( csv, "%s,%ld\n", "iparm_compress_reltol",     (long)iparm[IPARM_COMPRESS_RELTOL] );
    fprintf( csv, "%s,%ld\n", "iparm_compress_preselect",  (long)iparm[IPARM_COMPRESS_PRESELECT] );
    fprintf( csv, "%s,%ld\n", "iparm_compress_iluk",       (long)iparm[IPARM_COMPRESS_ILUK] );

    fprintf( csv, "%s,%ld\n", "iparm_mixed", (long)iparm[IPARM_MIXED] );
    fprintf( csv, "%s,%ld\n", "iparm_ftz",   (long)iparm[IPARM_FTZ] );

    fprintf( csv, "%s,%s\n",  "iparm_mpi_thread_level",  pastix_mpithreadmode_getstr(iparm[IPARM_MPI_THREAD_LEVEL]) );

    fprintf( csv, "%s,%ld\n", "iparm_modify_parameter", (long)iparm[IPARM_MODIFY_PARAMETER] );
    fprintf( csv, "%s,%s\n",  "iparm_start_task",        pastix_task_getstr(iparm[IPARM_START_TASK]) );
    fprintf( csv, "%s,%s\n",  "iparm_end_task",          pastix_task_getstr(iparm[IPARM_END_TASK]) );
    fprintf( csv, "%s,%s\n",  "iparm_float",             pastix_coeftype_getstr(iparm[IPARM_FLOAT]) );
    fprintf( csv, "%s,%ld\n", "iparm_mtx_type",         (long)iparm[IPARM_MTX_TYPE] );
    fprintf( csv, "%s,%ld\n", "iparm_dof_nbr",          (long)iparm[IPARM_DOF_NBR] );

    fprintf( csv, "%s,%e\n",  "dparm_fill_in",            dparm[DPARM_FILL_IN] );
    fprintf( csv, "%s,%e\n",  "dparm_epsilon_refinement", dparm[DPARM_EPSILON_REFINEMENT] );
    fprintf( csv, "%s,%e\n",  "dparm_relative_error",     dparm[DPARM_RELATIVE_ERROR] );
    fprintf( csv, "%s,%e\n",  "dparm_epsilon_magn_ctrl",  dparm[DPARM_EPSILON_MAGN_CTRL] );
    fprintf( csv, "%s,%e\n",  "dparm_order_time",         dparm[DPARM_ORDER_TIME] );
    fprintf( csv, "%s,%e\n",  "dparm_symbfact_time",      dparm[DPARM_SYMBFACT_TIME] );
    fprintf( csv, "%s,%e\n",  "dparm_reorder_time",       dparm[DPARM_REORDER_TIME] );
    fprintf( csv, "%s,%e\n",  "dparm_blend_time",         dparm[DPARM_BLEND_TIME] );
    fprintf( csv, "%s,%e\n",  "dparm_analyze_time",       dparm[DPARM_ANALYZE_TIME] );
    fprintf( csv, "%s,%e\n",  "dparm_pred_fact_time",     dparm[DPARM_PRED_FACT_TIME] );
    fprintf( csv, "%s,%e\n",  "dparm_fact_time",          dparm[DPARM_FACT_TIME] );
    fprintf( csv, "%s,%e\n",  "dparm_fact_flops",         dparm[DPARM_FACT_FLOPS] );
    fprintf( csv, "%s,%e\n",  "dparm_fact_thflops",       dparm[DPARM_FACT_THFLOPS] );
    fprintf( csv, "%s,%e\n",  "dparm_fact_rlflops",       dparm[DPARM_FACT_RLFLOPS] );
    fprintf( csv, "%s,%e\n",  "dparm_mem_fr",             dparm[DPARM_MEM_FR] );
    fprintf( csv, "%s,%e\n",  "dparm_mem_lr",             dparm[DPARM_MEM_LR] );
    fprintf( csv, "%s,%e\n",  "dparm_solv_time",          dparm[DPARM_SOLV_TIME] );
    fprintf( csv, "%s,%e\n",  "dparm_solv_flops",         dparm[DPARM_SOLV_FLOPS] );
    fprintf( csv, "%s,%e\n",  "dparm_solv_thflops",       dparm[DPARM_SOLV_THFLOPS] );
    fprintf( csv, "%s,%e\n",  "dparm_solv_rlflops",       dparm[DPARM_SOLV_RLFLOPS] );
    fprintf( csv, "%s,%e\n",  "dparm_refine_time",        dparm[DPARM_REFINE_TIME] );
    fprintf( csv, "%s,%e\n",  "dparm_a_norm",             dparm[DPARM_A_NORM] );
    fprintf( csv, "%s,%e\n",  "dparm_compress_tolerance", dparm[DPARM_COMPRESS_TOLERANCE] );
    fprintf( csv, "%s,%e\n",  "dparm_compress_min_ratio", dparm[DPARM_COMPRESS_MIN_RATIO] );
}
