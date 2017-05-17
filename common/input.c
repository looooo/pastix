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
    if(0 == strcmp("iparm_mc64", string))
    {
        return IPARM_MC64;
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
    if(0 == strcmp("iparm_compress_min_width", string))
    {
        return IPARM_COMPRESS_MIN_WIDTH;
    }
    if(0 == strcmp("iparm_compress_min_height", string))
    {
        return IPARM_COMPRESS_MIN_HEIGHT;
    }
    if(0 == strcmp("iparm_compress_when", string))
    {
        return IPARM_COMPRESS_WHEN;
    }
    if(0 == strcmp("iparm_compress_method", string))
    {
        return IPARM_COMPRESS_METHOD;
    }
    if(0 == strcmp("iparm_factorization", string))
    {
        return IPARM_FACTORIZATION;
    }
    if(0 == strcmp("iparm_scheduler", string))
    {
        return IPARM_SCHEDULER;
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
    if(0 == strcmp("iparm_gmres_im", string))
    {
        return IPARM_GMRES_IM;
    }
    if(0 == strcmp("iparm_free_cscuser", string))
    {
        return IPARM_FREE_CSCUSER;
    }
    if(0 == strcmp("iparm_autosplit_comm", string))
    {
        return IPARM_AUTOSPLIT_COMM;
    }
    if(0 == strcmp("iparm_gpu_nbr", string))
    {
        return IPARM_GPU_NBR;
    }
    if(0 == strcmp("iparm_produce_stats", string))
    {
        return IPARM_PRODUCE_STATS;
    }
    if(0 == strcmp("iparm_mtx_type", string))
    {
        return IPARM_MTX_TYPE;
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
    if(0 == strcmp("dparm_epsilon_refinement", string))
    {
        return DPARM_EPSILON_REFINEMENT;
    }
    if(0 == strcmp("dparm_epsilon_magn_ctrl", string))
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
    if(0 == strcmp("pastix_task_init", string))
    {
        return PastixTaskInit;
    }
    if( (0 == strcmp("pastixtaskordering", string)) ||
        (0 == strcmp("pastixtaskscotch",   string)) )
    {
        return PastixTaskOrdering;
    }
    if( (0 == strcmp("pastixtasksymbfact", string)) ||
        (0 == strcmp("pastixtaskfax",      string)) )
    {
        return PastixTaskSymbfact;
    }
    if( (0 == strcmp("pastixtaskanalyse", string)) ||
        (0 == strcmp("pastixtaskblend",   string)) )
    {
        return PastixTaskAnalyze;
    }
    if( (0 == strcmp("pastixtasknumfact", string)) ||
        (0 == strcmp("pastixtasksopalin", string)) )
    {
        return PastixTaskNumfact;
    }
    if( (0 == strcmp("pastixtasksolve",  string)) ||
        (0 == strcmp("pastixtaskupdown", string)) )
    {
        return PastixTaskSolve;
    }
    if( (0 == strcmp("pastixtaskrefine",     string)) ||
        (0 == strcmp("pastixtaskrefinement", string)) )
    {
        return PastixTaskRefine;
    }
    if(0 == strcmp("pastixtaskclean", string))
    {
        return PastixTaskClean;
    }
    if(0 == strcmp("pastixverbosenot", string))
    {
        return PastixVerboseNot;
    }
    if(0 == strcmp("pastixverboseno", string))
    {
        return PastixVerboseNo;
    }
    if(0 == strcmp("pastixverboseyes", string))
    {
        return PastixVerboseYes;
    }
    if(0 == strcmp("pastixiono", string))
    {
        return PastixIONo;
    }
    if(0 == strcmp("pastixioload", string))
    {
        return PastixIOLoad;
    }
    if(0 == strcmp("pastixiosave", string))
    {
        return PastixIOSave;
    }
    if(0 == strcmp("pastixioloadgraph", string))
    {
        return PastixIOLoadGraph;
    }
    if(0 == strcmp("pastixiosavegraph", string))
    {
        return PastixIOSaveGraph;
    }
    if(0 == strcmp("pastixioloadcsc", string))
    {
        return PastixIOLoadCSC;
    }
    if(0 == strcmp("pastixiosavecsc", string))
    {
        return PastixIOSaveCSC;
    }
    if(0 == strcmp("pastixrefinegmres", string))
    {
        return PastixRefineGMRES;
    }
    if(0 == strcmp("pastixrefinecg", string))
    {
        return PastixRefineCG;
    }
    if(0 == strcmp("pastixrefinesr", string))
    {
        return PastixRefineSR;
    }
    if(0 == strcmp("pastixrefinebicgstab", string))
    {
        return PastixRefineBiCGSTAB;
    }
    if(0 == strcmp("pastixorderscotch", string))
    {
        return PastixOrderScotch;
    }
    if(0 == strcmp("pastixordermetis", string))
    {
        return PastixOrderMetis;
    }
    if(0 == strcmp("pastixorderpersonal", string))
    {
        return PastixOrderPersonal;
    }
    if(0 == strcmp("pastixorderload", string))
    {
        return PastixOrderLoad;
    }
    if(0 == strcmp("pastixorderptscotch", string))
    {
        return PastixOrderPtscotch;
    }
    if(0 == strcmp("pastixorderparmetis", string))
    {
        return PastixOrderParMetis;
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
    if(0 == strcmp("pastixgeneral", string))
    {
        return PastixGeneral;
    }
    if(0 == strcmp("pastixhermitian", string))
    {
        return PastixHermitian;
    }
    if(0 == strcmp("pastixsymmetric", string))
    {
        return PastixSymmetric;
    }
    if(0 == strcmp("pastixschedsequential", string))
    {
        return PastixSchedSequential;
    }
    if(0 == strcmp("pastixschedstatic", string))
    {
        return PastixSchedStatic;
    }
    if(0 == strcmp("pastixscheddynamic", string))
    {
        return PastixSchedDynamic;
    }
    if(0 == strcmp("pastixschedparsec", string))
    {
        return PastixSchedParsec;
    }
    if(0 == strcmp("pastixschedstarpu", string))
    {
        return PastixSchedStarPU;
    }
    if(0 == strcmp("pastixcompressnever", string))
    {
        return PastixCompressNever;
    }
    if(0 == strcmp("pastixcompresswhenbegin", string))
    {
        return PastixCompressWhenBegin;
    }
    if(0 == strcmp("pastixcompresswhenend", string))
    {
        return PastixCompressWhenEnd;
    }
    if(0 == strcmp("pastixcompresswhenduring", string))
    {
        return PastixCompressWhenDuring;
    }
    if(0 == strcmp("pastixcompressmethodsvd", string))
    {
        return PastixCompressMethodSVD;
    }
    if(0 == strcmp("pastixcompressmethodrrqr", string))
    {
        return PastixCompressMethodRRQR;
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
