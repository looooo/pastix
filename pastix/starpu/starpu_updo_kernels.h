#ifndef STARPU_UPDO_KERNELS_H
#define STARPU_UPDO_KERNELS_H
#include "sopalin_define.h"
struct starpu_updo_trsm_data_ {
  pastix_int_t              cblknum;
  Sopalin_Data_t * sopalin_data;
  char             transpose;
  char             diag;
};
typedef struct starpu_updo_trsm_data_ starpu_updo_trsm_data_t;

struct starpu_updo_gemm_data_ {
  pastix_int_t              cblknum;
  pastix_int_t              bloknum;
  Sopalin_Data_t * sopalin_data;
  char             transpose;
};
typedef struct starpu_updo_gemm_data_ starpu_updo_gemm_data_t;

struct starpu_updo_diag_trsm_data_ {
  pastix_int_t              cblknum;
  Sopalin_Data_t * sopalin_data;
};
typedef struct starpu_updo_diag_trsm_data_ starpu_updo_diag_data_t;

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
#define updo_trsm_starpu_cpu API_CALL(updo_trsm_starpu_cpu)
void updo_trsm_starpu_cpu(void * buffers[], void * _args);
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
#define updo_down_gemm_starpu_cpu API_CALL(updo_down_gemm_starpu_cpu)
void updo_down_gemm_starpu_cpu(void * buffers[], void * _args);

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
#define updo_up_gemm_starpu_cpu API_CALL(updo_up_gemm_starpu_cpu)
void updo_up_gemm_starpu_cpu(void * buffers[], void * _args);

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
#define updo_diag_starpu_cpu API_CALL(updo_diag_starpu_cpu)
void updo_diag_starpu_cpu(void * buffers[], void * args);

#endif /* STARPU_UPDO_KERNELS_H */
