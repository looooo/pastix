#ifndef _OLD_API_H_
#define _OLD_API_H_

#define pastix_float_t void

/* Error numbers, need to conserve it MURGE compliant */
#define NO_ERR             PASTIX_SUCCESS
#define UNKNOWN_ERR        PASTIX_ERR_UNKNOWN
#define ALLOC_ERR          PASTIX_ERR_ALLOC
#define NOTIMPLEMENTED_ERR PASTIX_ERR_NOTIMPLEMENTED
#define OUTOFMEMORY_ERR    PASTIX_ERR_OUTOFMEMORY
#define THREAD_ERR         PASTIX_ERR_THREAD
#define INTERNAL_ERR       PASTIX_ERR_INTERNAL
#define BADPARAMETER_ERR   PASTIX_ERR_BADPARAMETER
#define FILE_ERR           PASTIX_ERR_FILE
#define INTEGER_TYPE_ERR   PASTIX_ERR_INTEGER_TYPE
#define IO_ERR             PASTIX_ERR_IO
#define MPI_ERR            PASTIX_ERR_MPI

/* Removed from PaStiX 6.0.0 */
#define ASSERT_ERR         -1
#define BAD_DEFINE_ERR     -1
#define FLOAT_TYPE_ERR     -1
#define MATRIX_ERR         -1
#define STEP_ORDER_ERR     -1

/* Former IPARM values */
/*
 * Backward compatibility
 */
enum IPARM_ACCESS_DEPRECATED {
    IPARM_DEFAULT_ORDERING      = IPARM_ORDERING_DEFAULT,
    IPARM_ORDERING_SWITCH_LEVEL = IPARM_SCOTCH_SWITCH_LEVEL,
    IPARM_ORDERING_CMIN         = IPARM_SCOTCH_CMIN,
    IPARM_ORDERING_CMAX         = IPARM_SCOTCH_CMAX,
    IPARM_ORDERING_FRAT         = IPARM_SCOTCH_FRAT,
    IPARM_AMALGAMATION_LEVEL    = IPARM_AMALGAMATION_LVLCBLK,
    IPARM_CUDA_NBR              = IPARM_GPU_NBR,
    IPARM_RHS_MAKING            = -1,
    IPARM_ONLY_RAFF             = IPARM_ONLY_REFINE,
    IPARM_MURGE_MAY_RAFF        = IPARM_MURGE_MAY_REFINE
};

/* Former DPARM values */
#define DPARM_RAFF_TIME DPARM_REFINE_TIME

/* Former API values */

/* _POS_ 1 */
enum API_TASK_OLD {
    API_TASK_SCOTCH     = API_TASK_ORDERING,
    API_TASK_FAX        = API_TASK_SYMBFACT,
    API_TASK_BLEND      = API_TASK_ANALYSE,
    API_TASK_SOPALIN    = API_TASK_NUMFACT,
    API_TASK_UPDOWN     = API_TASK_SOLVE,
    API_TASK_REFINEMENT = API_TASK_REFINE
};

/* _POS_ 4 */
#define API_FACT_LLT  PastixFactLLT
#define API_FACT_LDLT PastixFactLDLT
#define API_FACT_LU   PastixFactLU
#define API_FACT_LDLH PastixFactLDLH

/* Removed */
#define API_RHS_B -1
#define API_RHS_1 -1
#define API_RHS_I -1
#define API_RHS_0 -1

/* _POS_ 8 */
#define API_RAFF_GMRES    API_REFINE_GMRES
#define API_RAFF_GRAD     API_REFINE_GRAD
#define API_RAFF_PIVOT    API_REFINE_PIVOT
#define API_RAFF_BICGSTAB API_REFINE_BICGSTAB

/* _POS_ 61 */
#define API_REALSINGLE    PastixFloat
#define API_REALDOUBLE    PastixDouble
#define API_COMPLEXSINGLE PastixComplex32
#define API_COMPLEXDOUBLE PastixComplex64

/**
 * Some define for old pastix compatibility
 */
#define API_SYM_YES PastixSymmetric
#define API_SYM_HER PastixHermitian
#define API_SYM_NO  PastixGeneral

#endif /* _OLD_API_H_ */
