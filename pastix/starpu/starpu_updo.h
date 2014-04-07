/*
 * File: starpu_updo.h
 *
 * Updown step written using StarPU.
 *
 */
#ifndef STARPU_UPDO_H
#define STARPU_UPDO_H

#define starpu_register_sm2x API_CALL(starpu_register_sm2x)
int starpu_register_sm2x(Sopalin_Data_t       * sopalin_data,
                         starpu_data_handle_t * SM2X_handles);

#define starpu_submit_updown API_CALL(starpu_submit_updown)
int starpu_submit_updown(Sopalin_Data_t * sopalin_data,
                         starpu_data_handle_t * L_handles,
                         starpu_data_handle_t * U_handles,
                         starpu_data_handle_t * SM2X_handles,
			 struct starpu_task  ** tasktab,
                         int                  * sched_ctxs);
#endif /* STARPU_UPDO_H */
