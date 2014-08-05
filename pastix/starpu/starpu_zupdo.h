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
/*
 * File: starpu_updo.h
 *
 * Updown step written using StarPU.
 *
 */
#ifndef STARPU_UPDO_H
#define STARPU_UPDO_H

#define starpu_zregister_sm2x API_CALL(starpu_zregister_sm2x)
int starpu_zregister_sm2x(z_Sopalin_Data_t       * sopalin_data,
                         starpu_data_handle_t * SM2X_handles);

#define starpu_zsubmit_updown API_CALL(starpu_zsubmit_updown)
int starpu_zsubmit_updown(z_Sopalin_Data_t * sopalin_data,
                         starpu_data_handle_t * L_handles,
                         starpu_data_handle_t * U_handles,
                         starpu_data_handle_t * SM2X_handles,
			 struct starpu_task  ** tasktab,
                         int                  * sched_ctxs);
#endif /* STARPU_UPDO_H */
