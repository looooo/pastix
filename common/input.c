/**
 *
 * @file input.c
 *
 *  PaStiX API routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Théophile Terraz
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include <string.h>
#include "common.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_common
 *
 * iparm_to_int -
 *
 *******************************************************************************
 *
 * @param[in] string
 *          The iparm string to convert to int.
 *
 *******************************************************************************
 *
 * @return
 *      \retval value The returned integer.
 *
 *******************************************************************************/
int iparm_to_int( char * string)
{
    if(0 == strcmp("iparm_modify_parameter", string))
    {
        return IPARM_MODIFY_PARAMETER;
    }
    if(0 == strcmp("iparm_start_task", string))
    {
        return IPARM_START_TASK;
    }
    if(0 == strcmp("iparm_end_task", string))
    {
        return IPARM_END_TASK;
    }
    if(0 == strcmp("iparm_verbose", string))
    {
        return IPARM_VERBOSE;
    }
    if(0 == strcmp("iparm_dof_nbr", string))
    {
        return IPARM_DOF_NBR;
    }
    if(0 == strcmp("iparm_itermax", string))
    {
        return IPARM_ITERMAX;
    }
    if(0 == strcmp("iparm_matrix_verification", string))
    {
        return IPARM_MATRIX_VERIFICATION;
    }
    if(0 == strcmp("iparm_mc64", string))
    {
        return IPARM_MC64;
    }
    if(0 == strcmp("iparm_only_raff", string))
    {
        return IPARM_ONLY_RAFF;
    }
    if(0 == strcmp("iparm_cscd_correct", string))
    {
        return IPARM_CSCD_CORRECT;
    }
    if(0 == strcmp("iparm_tracefmt", string))
    {
        return IPARM_TRACEFMT;
    }
    if(0 == strcmp("iparm_ordering", string))
    {
        return IPARM_ORDERING;
    }
    if(0 == strcmp("iparm_ordering_default", string))
    {
        return IPARM_ORDERING_DEFAULT;
    }
    if(0 == strcmp("iparm_scotch_switch_level", string))
    {
        return IPARM_SCOTCH_SWITCH_LEVEL;
    }
    if(0 == strcmp("iparm_scotch_cmin", string))
    {
        return IPARM_SCOTCH_CMIN;
    }
    if(0 == strcmp("iparm_scotch_cmax", string))
    {
        return IPARM_SCOTCH_CMAX;
    }
    if(0 == strcmp("iparm_scotch_frat", string))
    {
        return IPARM_SCOTCH_FRAT;
    }
    if(0 == strcmp("iparm_metis_ctype", string))
    {
        return IPARM_METIS_CTYPE;
    }
    if(0 == strcmp("iparm_metis_rtype", string))
    {
        return IPARM_METIS_RTYPE;
    }
    if(0 == strcmp("iparm_metis_no2hop", string))
    {
        return IPARM_METIS_NO2HOP;
    }
    if(0 == strcmp("iparm_metis_nseps", string))
    {
        return IPARM_METIS_NSEPS;
    }
    if(0 == strcmp("iparm_metis_niter", string))
    {
        return IPARM_METIS_NITER;
    }
    if(0 == strcmp("iparm_metis_ufactor", string))
    {
        return IPARM_METIS_UFACTOR;
    }
    if(0 == strcmp("iparm_metis_compress", string))
    {
        return IPARM_METIS_COMPRESS;
    }
    if(0 == strcmp("iparm_metis_ccorder", string))
    {
        return IPARM_METIS_CCORDER;
    }
    if(0 == strcmp("iparm_metis_pfactor", string))
    {
        return IPARM_METIS_PFACTOR;
    }
    if(0 == strcmp("iparm_metis_seed", string))
    {
        return IPARM_METIS_SEED;
    }
    if(0 == strcmp("iparm_metis_dbglvl", string))
    {
        return IPARM_METIS_DBGLVL;
    }
    if(0 == strcmp("iparm_sf_kass", string))
    {
        return IPARM_SF_KASS;
    }
    if(0 == strcmp("iparm_amalgamation_lvlcblk", string))
    {
        return IPARM_AMALGAMATION_LVLCBLK;
    }
    if(0 == strcmp("iparm_amalgamation_lvlblas", string))
    {
        return IPARM_AMALGAMATION_LVLBLAS;
    }
    if(0 == strcmp("iparm_reordering_split", string))
    {
        return IPARM_REORDERING_SPLIT;
    }
    if(0 == strcmp("iparm_reordering_stop", string))
    {
        return IPARM_REORDERING_STOP;
    }
    if(0 == strcmp("iparm_baseval", string))
    {
        return IPARM_BASEVAL;
    }
    if(0 == strcmp("iparm_min_blocksize", string))
    {
        return IPARM_MIN_BLOCKSIZE;
    }
    if(0 == strcmp("iparm_max_blocksize", string))
    {
        return IPARM_MAX_BLOCKSIZE;
    }
    if(0 == strcmp("iparm_compress_size", string))
    {
        return IPARM_COMPRESS_SIZE;
    }
    if(0 == strcmp("iparm_compress_width", string))
    {
        return IPARM_COMPRESS_WIDTH;
    }
    if(0 == strcmp("iparm_compress_when", string))
    {
        return IPARM_COMPRESS_WHEN;
    }
    if(0 == strcmp("iparm_compress_method", string))
    {
        return IPARM_COMPRESS_METHOD;
    }
    if(0 == strcmp("iparm_schur", string))
    {
        return IPARM_SCHUR;
    }
    if(0 == strcmp("iparm_isolate_zeros", string))
    {
        return IPARM_ISOLATE_ZEROS;
    }
    if(0 == strcmp("iparm_rhsd_check", string))
    {
        return IPARM_RHSD_CHECK;
    }
    if(0 == strcmp("iparm_factorization", string))
    {
        return IPARM_FACTORIZATION;
    }
    if(0 == strcmp("iparm_scheduler", string))
    {
        return IPARM_SCHEDULER;
    }
    if(0 == strcmp("iparm_cpu_by_node", string))
    {
        return IPARM_CPU_BY_NODE;
    }
    if(0 == strcmp("iparm_bindthrd", string))
    {
        return IPARM_BINDTHRD;
    }
    if(0 == strcmp("iparm_thread_nbr", string))
    {
        return IPARM_THREAD_NBR;
    }
    if(0 == strcmp("iparm_distribution_level", string))
    {
        return IPARM_DISTRIBUTION_LEVEL;
    }
    if(0 == strcmp("iparm_level_of_fill", string))
    {
        return IPARM_LEVEL_OF_FILL;
    }
    if(0 == strcmp("iparm_io_strategy", string))
    {
        return IPARM_IO_STRATEGY;
    }
    if(0 == strcmp("iparm_rhs_making", string))
    {
        return IPARM_RHS_MAKING;
    }
    if(0 == strcmp("iparm_refinement", string))
    {
        return IPARM_REFINEMENT;
    }
    if(0 == strcmp("iparm_incomplete", string))
    {
        return IPARM_INCOMPLETE;
    }
    if(0 == strcmp("iparm_abs", string))
    {
        return IPARM_ABS;
    }
    if(0 == strcmp("iparm_esp", string))
    {
        return IPARM_ESP;
    }
    if(0 == strcmp("iparm_gmres_im", string))
    {
        return IPARM_GMRES_IM;
    }
    if(0 == strcmp("iparm_free_cscuser", string))
    {
        return IPARM_FREE_CSCUSER;
    }
    if(0 == strcmp("iparm_free_cscpastix", string))
    {
        return IPARM_FREE_CSCPASTIX;
    }
    if(0 == strcmp("iparm_ooc_limit", string))
    {
        return IPARM_OOC_LIMIT;
    }
    if(0 == strcmp("iparm_ooc_thread", string))
    {
        return IPARM_OOC_THREAD;
    }
    if(0 == strcmp("iparm_ooc_id", string))
    {
        return IPARM_OOC_ID;
    }
    if(0 == strcmp("iparm_nb_smp_node_used", string))
    {
        return IPARM_NB_SMP_NODE_USED;
    }
    if(0 == strcmp("iparm_thread_comm_mode", string))
    {
        return IPARM_THREAD_COMM_MODE;
    }
    if(0 == strcmp("iparm_nb_thread_comm", string))
    {
        return IPARM_NB_THREAD_COMM;
    }
    if(0 == strcmp("iparm_fill_matrix", string))
    {
        return IPARM_FILL_MATRIX;
    }
    if(0 == strcmp("iparm_esp_threshold", string))
    {
        return IPARM_ESP_THRESHOLD;
    }
    if(0 == strcmp("iparm_dof_cost", string))
    {
        return IPARM_DOF_COST;
    }
    if(0 == strcmp("iparm_murge_refinement", string))
    {
        return IPARM_MURGE_REFINEMENT;
    }
    if(0 == strcmp("iparm_starpu", string))
    {
        return IPARM_STARPU;
    }
    if(0 == strcmp("iparm_autosplit_comm", string))
    {
        return IPARM_AUTOSPLIT_COMM;
    }
    if(0 == strcmp("iparm_cuda_nbr", string))
    {
        return IPARM_CUDA_NBR;
    }
    if(0 == strcmp("iparm_transpose_solve", string))
    {
        return IPARM_TRANSPOSE_SOLVE;
    }
    if(0 == strcmp("iparm_starpu_ctx_depth", string))
    {
        return IPARM_STARPU_CTX_DEPTH;
    }
    if(0 == strcmp("iparm_starpu_ctx_nbr", string))
    {
        return IPARM_STARPU_CTX_NBR;
    }
    if(0 == strcmp("iparm_produce_stats", string))
    {
        return IPARM_PRODUCE_STATS;
    }
    if(0 == strcmp("iparm_gpu_criterium", string))
    {
        return IPARM_GPU_CRITERIUM;
    }
    if(0 == strcmp("iparm_murge_may_refine", string))
    {
        return IPARM_MURGE_MAY_REFINE;
    }
    return -1;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_common
 *
 * dparm_to_int -
 *
 *******************************************************************************
 *
 * @param[in] string
 *          The dparm string to convert to int.
 *
 *******************************************************************************
 *
 * @return
 *      \retval value The returned integer.
 *
 *******************************************************************************/
int dparm_to_int( char * string)
{
    if(0 == strcmp("DPARM_EPSILON_REFINEMENT", string))
    {
        return DPARM_EPSILON_REFINEMENT;
    }
    if(0 == strcmp("DPARM_EPSILON_MAGN_CTRL", string))
    {
        return DPARM_EPSILON_MAGN_CTRL;
    }
    if(0 == strcmp("dparm_compress_tolerance", string))
    {
        return DPARM_COMPRESS_TOLERANCE;
    }
    return -1;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_common
 *
 * api_to_int -
 *
 *******************************************************************************
 *
 * @param[in] string
 *          The api string to convert to int.
 *
 *******************************************************************************
 *
 * @return
 *      \retval value The returned integer.
 *
 *******************************************************************************/
int api_to_int( char * string )
{
    if(0 == strcmp("api_task_init", string))
    {
        return API_TASK_INIT;
    }
    if(0 == strcmp("api_task_ordering", string))
    {
        return API_TASK_ORDERING;
    }
    if(0 == strcmp("api_task_symbfact", string))
    {
        return API_TASK_SYMBFACT;
    }
    if(0 == strcmp("api_task_analyse", string))
    {
        return API_TASK_ANALYSE;
    }
    if(0 == strcmp("api_task_numfact", string))
    {
        return API_TASK_NUMFACT;
    }
    if(0 == strcmp("api_task_solve", string))
    {
        return API_TASK_SOLVE;
    }
    if(0 == strcmp("api_task_refine", string))
    {
        return API_TASK_REFINE;
    }
    if(0 == strcmp("api_task_clean", string))
    {
        return API_TASK_CLEAN;
    }
    if(0 == strcmp("api_task_scotch", string))
    {
        return API_TASK_SCOTCH;
    }
    if(0 == strcmp("api_task_fax", string))
    {
        return API_TASK_FAX;
    }
    if(0 == strcmp("api_task_blend", string))
    {
        return API_TASK_BLEND;
    }
    if(0 == strcmp("api_task_sopalin", string))
    {
        return API_TASK_SOPALIN;
    }
    if(0 == strcmp("api_task_updown", string))
    {
        return API_TASK_UPDOWN;
    }
    if(0 == strcmp("api_task_refinement", string))
    {
        return API_TASK_REFINEMENT;
    }
    if(0 == strcmp("api_verbose_not", string))
    {
        return API_VERBOSE_NOT;
    }
    if(0 == strcmp("api_verbose_no", string))
    {
        return API_VERBOSE_NO;
    }
    if(0 == strcmp("api_verbose_yes", string))
    {
        return API_VERBOSE_YES;
    }
    if(0 == strcmp("api_verbose_chatterbox", string))
    {
        return API_VERBOSE_CHATTERBOX;
    }
    if(0 == strcmp("api_verbose_unbearable", string))
    {
        return API_VERBOSE_UNBEARABLE;
    }
    if(0 == strcmp("api_io_no", string))
    {
        return API_IO_NO;
    }
    if(0 == strcmp("api_io_load", string))
    {
        return API_IO_LOAD;
    }
    if(0 == strcmp("api_io_save", string))
    {
        return API_IO_SAVE;
    }
    if(0 == strcmp("api_io_load_graph", string))
    {
        return API_IO_LOAD_GRAPH;
    }
    if(0 == strcmp("api_io_save_graph", string))
    {
        return API_IO_SAVE_GRAPH;
    }
    if(0 == strcmp("api_io_load_csc", string))
    {
        return API_IO_LOAD_CSC;
    }
    if(0 == strcmp("api_io_save_csc", string))
    {
        return API_IO_SAVE_CSC;
    }
    if(0 == strcmp("api_rhs_b", string))
    {
        return API_RHS_B;
    }
    if(0 == strcmp("api_rhs_1", string))
    {
        return API_RHS_1;
    }
    if(0 == strcmp("api_rhs_i", string))
    {
        return API_RHS_I;
    }
    if(0 == strcmp("api_rhs_0", string))
    {
        return API_RHS_0;
    }
    if(0 == strcmp("api_raf_gmres", string))
    {
        return API_RAF_GMRES;
    }
    if(0 == strcmp("api_raf_grad", string))
    {
        return API_RAF_GRAD;
    }
    if(0 == strcmp("api_raf_pivot", string))
    {
        return API_RAF_PIVOT;
    }
    if(0 == strcmp("api_raf_bicgstab", string))
    {
        return API_RAF_BICGSTAB;
    }
    if(0 == strcmp("api_csc_preserve", string))
    {
        return API_CSC_PRESERVE;
    }
    if(0 == strcmp("api_csc_free", string))
    {
        return API_CSC_FREE;
    }
    if(0 == strcmp("api_thread_multiple", string))
    {
        return API_THREAD_MULTIPLE;
    }
    if(0 == strcmp("api_thread_funneled", string))
    {
        return API_THREAD_FUNNELED;
    }
    if(0 == strcmp("api_thread_comm_one", string))
    {
        return API_THREAD_COMM_ONE;
    }
    if(0 == strcmp("api_thread_comm_defined", string))
    {
        return API_THREAD_COMM_DEFINED;
    }
    if(0 == strcmp("api_thread_comm_nbproc", string))
    {
        return API_THREAD_COMM_NBPROC;
    }
    if(0 == strcmp("api_bind_no", string))
    {
        return API_BIND_NO;
    }
    if(0 == strcmp("api_bind_auto", string))
    {
        return API_BIND_AUTO;
    }
    if(0 == strcmp("api_bind_tab", string))
    {
        return API_BIND_TAB;
    }
    if(0 == strcmp("api_no", string))
    {
        return API_NO;
    }
    if(0 == strcmp("api_yes", string))
    {
        return API_YES;
    }
    if(0 == strcmp("api_trace_paje", string))
    {
        return API_TRACE_PAJE;
    }
    if(0 == strcmp("api_trace_humread", string))
    {
        return API_TRACE_HUMREAD;
    }
    if(0 == strcmp("api_trace_unformated", string))
    {
        return API_TRACE_UNFORMATED;
    }
    if(0 == strcmp("api_order_scotch", string))
    {
        return API_ORDER_SCOTCH;
    }
    if(0 == strcmp("api_order_metis", string))
    {
        return API_ORDER_METIS;
    }
    if(0 == strcmp("api_order_personal", string))
    {
        return API_ORDER_PERSONAL;
    }
    if(0 == strcmp("api_order_load", string))
    {
        return API_ORDER_LOAD;
    }
    if(0 == strcmp("api_order_ptscotch", string))
    {
        return API_ORDER_PTSCOTCH;
    }
    if(0 == strcmp("api_gpu_criterion_updates", string))
    {
        return API_GPU_CRITERION_UPDATES;
    }
    if(0 == strcmp("api_gpu_criterion_cblksize", string))
    {
        return API_GPU_CRITERION_CBLKSIZE;
    }
    if(0 == strcmp("api_gpu_criterion_flops", string))
    {
        return API_GPU_CRITERION_FLOPS;
    }
    if(0 == strcmp("api_gpu_criterion_priority", string))
    {
        return API_GPU_CRITERION_PRIORITY;
    }
    if(0 == strcmp("pastixfactllt", string))
    {
        return PastixFactLLT;
    }
    if(0 == strcmp("pastixfactldlt", string))
    {
        return PastixFactLDLT;
    }
    if(0 == strcmp("pastixfactlu", string))
    {
        return PastixFactLU;
    }
    if(0 == strcmp("pastixfactldlh", string))
    {
        return PastixFactLDLH;
    }
    if(atoi(string) == 0)
    {
        if(0 == strcmp("0", string))
        {
            return 0;
        }
        else
        {
            return -1;
        }
    }
    else
    {
        return atoi(string);
    }
}
