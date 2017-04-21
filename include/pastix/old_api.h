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

#define ASSERT_ERR         -1
#define BAD_DEFINE_ERR     -1
#define FLOAT_TYPE_ERR     -1
#define MATRIX_ERR         -1
#define STEP_ORDER_ERR     -1

/* Former IPARM values */
#define IPARM_ONLY_RAFF      IPARM_ONLY_REFINE
#define IPARM_MURGE_MAY_RAFF IPARM_MURGE_MAY_REFINE

/* Former DPARM values */
#define DPARM_RAFF_TIME DPARM_REFINE_TIME

/* Former API values */

/* _POS_ 1 */
#define API_TASK_RAFF API_TASK_REFINE

/* _POS_ 4 */
#define API_FACT_LLT  PastixFactLLT
#define API_FACT_LDLT PastixFactLDLT
#define API_FACT_LU   PastixFactLU
#define API_FACT_LDLH PastixFactLDLH

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
