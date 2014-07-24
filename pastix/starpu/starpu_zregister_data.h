/**
 * @file starpu_zregister_data.h
 *
 * @author Xavier Lacoste
 * @precisions normal z -> s d c
 */

#include "common.h"
#include "starpu_zdefines.h"
#include "z_solver.h"


void starpu_zfanin_init_cpu_func(void *descr[], void *cl_arg);

int
starpu_zregister_fanin(z_SolverMatrix            * solvmtx,
                       starpu_data_handle_t  *** Lfanin_handle,
                       starpu_data_handle_t  *** Ufanin_handle);

int
starpu_zregister_halo(z_SolverMatrix          * datacode,
                      starpu_data_handle_t ** Lhalo_handle,
                      starpu_data_handle_t ** Uhalo_handle);

int
starpu_zregister_cblk( z_SolverMatrix          * datacode,
                       starpu_data_handle_t ** L_handle,
                       starpu_data_handle_t ** U_handle );

int
starpu_zregister_blocktab( Sopalin_Data_t        * sopalin_data,
                           starpu_data_handle_t ** blocktab_handles,
                           int                  ** blocktab);

int
starpu_zregister_work( z_SolverMatrix * datacode,
                       starpu_data_handle_t * WORK_handle,
                       pastix_int_t WORK_size );

int
starpu_zregister_data( Sopalin_Data_t         * sopalin_data,
                       starpu_data_handle_t  ** L_handle,
                       starpu_data_handle_t  ** U_handle,
                       starpu_data_handle_t  ** Lhalo_handle,
                       starpu_data_handle_t  ** Uhalo_handle,
                       starpu_data_handle_t *** Lfanin_handle,
                       starpu_data_handle_t *** Ufanin_handle,
                       starpu_data_handle_t  ** blocktab_handles,
                       int                   ** blocktab,
                       starpu_data_handle_t   * WORK_handle,
                       pastix_int_t             WORK_size );

int
starpu_zunregister_fanin(z_SolverMatrix            * solvmtx,
                         starpu_data_handle_t  *** Lfanin_handle,
                         starpu_data_handle_t  *** Ufanin_handle);

int
starpu_zunregister_halo(z_SolverMatrix          * datacode,
                        starpu_data_handle_t ** Lhalo_handle,
                        starpu_data_handle_t ** Uhalo_handle);

int
starpu_zunregister_cblk( z_SolverMatrix          * datacode,
                         starpu_data_handle_t ** L_handle,
                         starpu_data_handle_t ** U_handle );

int
starpu_zunregister_blocktab( z_SolverMatrix          * datacode,
                             starpu_data_handle_t ** blocktab_handles,
                             int                  ** blocktab);

int
starpu_zunregister_work( z_SolverMatrix * datacode,
                         starpu_data_handle_t * WORK_handle );

int
starpu_zunregister_data( Sopalin_Data_t         * sopalin_data,
                         starpu_data_handle_t  ** L_handle,
                         starpu_data_handle_t  ** U_handle,
                         starpu_data_handle_t  ** Lhalo_handle,
                         starpu_data_handle_t  ** Uhalo_handle,
                         starpu_data_handle_t *** Lfanin_handle,
                         starpu_data_handle_t *** Ufanin_handle,
                         starpu_data_handle_t  ** blocktab_handles,
                         int                   ** blocktab,
                         starpu_data_handle_t   * WORK_handle);
