/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#ifndef Z_OOC_H
#define Z_OOC_H

/* Return values */
#define EXIT_FAILURE_CBLK_NOT_NULL      3
#define EXIT_FAILURE_SAVING_NULL_BUFFER 4
#define EXIT_FAILURE_OUT_OF_MEMORY      5
#define EXIT_FAILURE_FILE_OPENING       6
#define EXIT_FAILURE_FILE_TRUNCATED     7
#define EXIT_SUCCESS_HACK               8
#define EXIT_SUCCESS_ALL_LOADED         9
#define EXIT_FAILURE_CBLK_USED          10

/* OOC Step */
#define OOCSTEP_COEFINIT                1
#define OOCSTEP_SOPALIN                 2
#define OOCSTEP_DOWN                    3
#define OOCSTEP_DIAG                    4
#define OOCSTEP_UP                      5


#ifdef OOC

#define OOC_RECEIVING z_ooc_receiving(sopalin_data)
#define OOC_RECEIVED z_ooc_received(sopalin_data)
#define OOC_THREAD_NBR sopalin_data->sopar->iparm[IPARM_OOC_THREAD]

/* sets values for the global ooc structure 
 * sopalin_data : z_Sopalin_Data_t global structure
 * limit        : memory limit set by use
 */
void *z_ooc_thread(void * arg);

/* Init / Clean */
int z_ooc_init          (z_Sopalin_Data_t * sopalin_data, pastix_int_t limit);
int z_ooc_exit          (z_Sopalin_Data_t * sopalin_data);

/* Step */
int z_ooc_stop_thread   (z_Sopalin_Data_t * sopalin_data);
int z_ooc_freeze        (z_Sopalin_Data_t * sopalin_data);
int z_ooc_defreeze      (z_Sopalin_Data_t * sopalin_data);
int z_ooc_set_step      (z_Sopalin_Data_t * sopalin_data, int step);

/* Cblk */
int z_ooc_wait_for_cblk (z_Sopalin_Data_t * sopalin_data, pastix_int_t cblk, int me);
int z_ooc_hack_load     (z_Sopalin_Data_t * sopalin_data, pastix_int_t cblk, int me);
int z_ooc_save_coef     (z_Sopalin_Data_t * sopalin_data, pastix_int_t task, pastix_int_t cblk, int me);

void z_ooc_receiving    (z_Sopalin_Data_t * sopalin_data);
void z_ooc_received     (z_Sopalin_Data_t * sopalin_data);
void z_ooc_wait_task    (z_Sopalin_Data_t * sopalin_data, pastix_int_t task, int me);
#else /* OOC */

#define OOC_RECEIVING 
#define OOC_RECEIVED 
#define OOC_THREAD_NBR 0
#define z_ooc_thread     NULL

#define z_ooc_init(sopalin_data, limit)
#define z_ooc_exit(sopalin_data)

#define z_ooc_stop_thread(sopalin_data)
#define z_ooc_freeze(sopalin_data)
#define z_ooc_defreeze(sopalin_data)
#define z_ooc_set_step(sopalin_data, step)

#define z_ooc_wait_for_cblk(sopalin_data, cblk, me)
#define z_ooc_save_coef(sopalin_data, task, cblk, me)
#define z_ooc_hack_load(sopalin_data, cblknum, me)

#define z_ooc_wait_task(sopalin_data, task, me) 

#endif /* OOC */

#ifdef OOC_FTGT
/* Ftgt */
int z_ooc_wait_for_ftgt (z_Sopalin_Data_t * sopalin_data, pastix_int_t ftgtnum, int me);
int z_ooc_reset_ftgt    (z_Sopalin_Data_t * sopalin_data, pastix_int_t ftgtnum, int me);
int z_ooc_save_ftgt     (z_Sopalin_Data_t * sopalin_data, pastix_int_t tasknum, pastix_int_t ftgtnum, int me);

#else

#define z_ooc_wait_for_ftgt(sopalin_data, ftgtnum, me)
#define z_ooc_reset_ftgt(sopalin_data, ftgtnum, me)
#define z_ooc_save_ftgt(sopalin_data, tasknum, ftgtnum, me)
#endif

#endif /* Z_OOC_H */
