/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.2.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#ifdef PASTIX_WITH_MAGMABLAS
#  include <magmablas.h>
#endif /* PASTIX_WITH_MAGMABLAS */


#include <starpu.h>
#include "common.h"
#include "z_sopalin3d.h"
#include "sopalin_acces.h"
#include "starpu_zkernels.h"
#include "starpu_zupdo_kernels.h"
#include "z_compute_trsm.h"
#include <inttypes.h>



/*
 * Function: updo_trsm_starpu_common
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     1            - L column block
 *     2            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
static inline
void updo_trsm_starpu_common(void * buffers[], void * _args, int arch)
{
  z_Sopalin_Data_t          * sopalin_data;
  z_SolverMatrix            * datacode;
  pastix_complex64_t            * L    = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
  pastix_complex64_t            * RHS  = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
  pastix_int_t                stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  pastix_int_t                rhsnbr       = STARPU_MATRIX_GET_NY(buffers[1]);
  pastix_int_t                rhssze       = STARPU_MATRIX_GET_LD(buffers[1]);
  pastix_int_t                cblknum;
  char                      transpose;
  char                      diag;
  pastix_int_t                colnbr;
  pastix_complex64_t              fun          = 1.0;
  starpu_codelet_unpack_args(_args, &sopalin_data, &cblknum,
                             &transpose, &diag);
  datacode = sopalin_data->datacode;
  colnbr = CBLK_COLNBR(cblknum);
  ASSERTDBG(UPDOWN_SM2XNBR == rhsnbr, MOD_SOPALIN);
  ASSERTDBG(UPDOWN_SM2XSZE == rhssze, MOD_SOPALIN);
  switch(arch) {
  case ARCH_CPU:
    SOPALIN_TRSM("L","L", &transpose, &diag,colnbr,rhsnbr,fun,L,stride,RHS,rhssze);
    break;
  case ARCH_CUDA:
  default:
    errorPrint("Unknown Architecture");
    assert(0);
    break;
  }
}



/*
 * Function: updo_trsm_starpu_cpu
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     1            - L column block
 *     2            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
void updo_trsm_starpu_cpu(void * buffers[], void * _args)
{
  updo_trsm_starpu_common(buffers, _args, ARCH_CPU);
}


/*
 * Function: updo_down_gemm_starpu_common
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     1            - L column block
 *     2            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
static inline
void updo_down_gemm_starpu_common(void * buffers[], void * _args, int arch)
{
  z_Sopalin_Data_t          * sopalin_data;
  z_SolverMatrix            * datacode;
  pastix_complex64_t            * L            = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
  pastix_complex64_t            * RHS          = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
  pastix_complex64_t            * RHS2         = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[2]);
  pastix_int_t                stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  pastix_int_t                rhsnbr       = STARPU_MATRIX_GET_NY(buffers[1]);
  pastix_int_t                rhssze       = STARPU_MATRIX_GET_LD(buffers[1]);
  pastix_int_t                cblknum;
  pastix_int_t                bloknum;
  char                      transpose;
  pastix_int_t                fcblknum;
  pastix_int_t                colnbr;
  pastix_int_t                rownbr;
  pastix_complex64_t              fun          = 1.0;
  pastix_complex64_t            * ga;
  pastix_complex64_t            * gc;

  starpu_codelet_unpack_args(_args, &sopalin_data, &cblknum, &bloknum,
                             &transpose);
  datacode = sopalin_data->datacode;
  fcblknum = SYMB_CBLKNUM(bloknum);
  colnbr = CBLK_COLNBR(cblknum);
  rownbr = BLOK_ROWNBR(bloknum);
  ga = L + SOLV_COEFIND(bloknum);
  gc = RHS2 +
    SYMB_FROWNUM(bloknum) - SYMB_FCOLNUM(fcblknum);

  ASSERTDBG(UPDOWN_SM2XNBR == rhsnbr, MOD_SOPALIN);
  ASSERTDBG(UPDOWN_SM2XSZE == rhssze, MOD_SOPALIN);

  switch(arch) {
  case ARCH_CPU:
    SOPALIN_GEMM(&transpose,"N",rownbr,rhsnbr,colnbr,-fun,ga,stride,
                 RHS,rhssze,fun,gc, UPDOWN_SM2XSZE);
    break;
  case ARCH_CUDA:
  default:
    errorPrint("Unknown Architecture");
    assert(0);
    break;
  }
}



/*
 * Function: updo_down_gemm_starpu_cpu
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     1            - L column block
 *     2            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
void updo_down_gemm_starpu_cpu(void * buffers[], void * _args)
{
  updo_down_gemm_starpu_common(buffers, _args, ARCH_CPU);
}


/*
 * Function: updo_up_gemm_starpu_common
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     1            - L column block
 *     2            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
static inline
void updo_up_gemm_starpu_common(void * buffers[], void * _args, int arch)
{
  z_Sopalin_Data_t          * sopalin_data;
  z_SolverMatrix            * datacode;
  pastix_complex64_t            * L            = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
  pastix_complex64_t            * RHS          = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
  pastix_complex64_t            * RHS2         = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[2]);
  pastix_int_t                stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  pastix_int_t                rhsnbr       = STARPU_MATRIX_GET_NY(buffers[1]);
  pastix_int_t                rhssze       = STARPU_MATRIX_GET_LD(buffers[1]);
  pastix_int_t                cblknum;
  pastix_int_t                bloknum;
  char                      transpose;
  pastix_int_t                fcblknum;
  pastix_int_t                colnbr;
  pastix_int_t                rownbr;
  pastix_complex64_t              fun          = 1.0;
  pastix_complex64_t            * ga;
  pastix_complex64_t            * gc;

  starpu_codelet_unpack_args(_args, &sopalin_data, &cblknum, &bloknum,
                             &transpose);
  datacode = sopalin_data->datacode;
  fcblknum = SYMB_CBLKNUM(bloknum);
  colnbr = CBLK_COLNBR(cblknum);
  rownbr = BLOK_ROWNBR(bloknum);
  ga = L + SOLV_COEFIND(bloknum);
  gc = RHS2 +
    SYMB_FROWNUM(bloknum) - SYMB_FCOLNUM(fcblknum);

  ASSERTDBG(UPDOWN_SM2XNBR == rhsnbr, MOD_SOPALIN);
  ASSERTDBG(UPDOWN_SM2XSZE == rhssze, MOD_SOPALIN);

  switch(arch) {
  case ARCH_CPU:
    SOPALIN_GEMM(&transpose,"N",colnbr,rhsnbr,rownbr,-fun,ga,stride,
                 gc,rhssze,fun,RHS, UPDOWN_SM2XSZE);
    break;
  case ARCH_CUDA:
  default:
    errorPrint("Unknown Architecture");
    assert(0);
    break;
  }
}



/*
 * Function: updo_up_gemm_starpu_cpu
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     1            - L column block
 *     2            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
void updo_up_gemm_starpu_cpu(void * buffers[], void * _args)
{
  updo_up_gemm_starpu_common(buffers, _args, ARCH_CPU);
}

/*
 * Function: updo_diag_starpu_common
 *
 * Divide the right-hand-side(s) by the diagonal
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     0            - L column block
 *     1            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 *   arch       - Type of architecture : ARCH_CPU | ARCH_CUDA
 */
static inline
void updo_diag_starpu_common(void * buffers[], void * _args, int arch)
{
  z_Sopalin_Data_t          * sopalin_data;
  z_SolverMatrix            * datacode;
  pastix_complex64_t            * L            = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
  pastix_complex64_t            * RHS          = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
  pastix_int_t                stride       = STARPU_MATRIX_GET_LD(buffers[0]);
  pastix_int_t                rhsnbr       = STARPU_MATRIX_GET_NY(buffers[1]);
  pastix_int_t                rhssze       = STARPU_MATRIX_GET_LD(buffers[1]);
  pastix_int_t                cblknum;
  pastix_int_t                colnbr;

  starpu_codelet_unpack_args(_args, &sopalin_data, &cblknum);
  datacode = sopalin_data->datacode;
  colnbr = CBLK_COLNBR(cblknum);

  ASSERTDBG(UPDOWN_SM2XNBR == rhsnbr, MOD_SOPALIN);
  ASSERTDBG(UPDOWN_SM2XSZE == rhssze, MOD_SOPALIN);

  switch(arch) {
  case ARCH_CPU:
  {
    pastix_int_t i, j;
    pastix_complex64_t * myRHS = RHS;
    for (j = 0; j < rhsnbr; j++)
      {
        for (i = 0; i < colnbr; i++)
          {
            myRHS[i] /= L[i*(stride+1)];
          }
        myRHS += rhssze;
      }
    break;
  }
  case ARCH_CUDA:
  default:
    errorPrint("Unknown Architecture");
    assert(0);
    break;
  }
}

/*
 * Function: updo_diag_starpu_cpu
 *
 * Divide the right-hand-side(s) by the diagonal.
 *
 * CPU interface.
 *
 * Parameters:
 *   buffers    - Data handlers :
 *     0            - L column block
 *     1            - Right-hand-side block facing the column block.
 *   _args      - Codelet arguments:
 *     sopalin_data - global PaStiX internal data.
 *     cblknum      - Current column block index.
 */
void updo_diag_starpu_cpu(void * buffers[], void * args)
{
  updo_diag_starpu_common(buffers, args, ARCH_CPU);
}
