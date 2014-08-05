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
#ifndef Z_STARPU_ZSUBMIT_TASKS_H
#define Z_STARPU_ZSUBMIT_TASKS_H

#define starpu_zloop_data_  API_CALL(starpu_zloop_data_)
#define starpu_zloop_data_t API_CALL(starpu_zloop_data_t)
typedef struct starpu_zloop_data_ starpu_zloop_data_t;

#define starpu_zsubmit_one_trf API_CALL(starpu_zsubmit_one_trf)
int starpu_zsubmit_one_trf (pastix_int_t itertask, z_Sopalin_Data_t * sopalin_data);
#define starpu_zsubmit_bunch_of_gemm API_CALL(starpu_zsubmit_bunch_of_gemm)
int starpu_zsubmit_bunch_of_gemm (pastix_int_t itertask, z_Sopalin_Data_t * sopalin_data);

#define starpu_zsubmit_tasks API_CALL(starpu_zsubmit_tasks)
int starpu_zsubmit_tasks(z_Sopalin_Data_t   * sopalin_data);

#endif /* Z_STARPU_ZSUBMIT_TASKS_H */
