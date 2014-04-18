/**
 * @file starpu_zregister_data.h
 *
 * @author Xavier Lacoste
 * @precisions normal z -> s d c
 */

#include "common.h"
#include "starpu_defines.h"
#include "solver.h"


void starpu_zfanin_init_cpu_func(void *descr[], void *cl_arg);

int
starpu_zregister_fanin(SolverMatrix            * solvmtx,
                       starpu_data_handle_t  *** Lfanin_handle,
                       starpu_data_handle_t  *** Ufanin_handle);

int
starpu_zregister_halo(SolverMatrix          * datacode,
                      starpu_data_handle_t ** Lhalo_handle,
                      starpu_data_handle_t ** Uhalo_handle);

int
starpu_zregister_cblk( SolverMatrix          * datacode,
                       starpu_data_handle_t ** L_handle,
                       starpu_data_handle_t ** U_handle );

int
starpu_zregister_data( Sopalin_Data_t * sopalin_data,
                       starpu_data_handle_t ** L_handle,
                       starpu_data_handle_t ** U_handle,
                       starpu_data_handle_t ** Lhalo_handle,
                       starpu_data_handle_t ** Uhalo_handle,
                       starpu_data_handle_t *** Lfanin_handle,
                       starpu_data_handle_t *** Ufanin_handle);

int
starpu_zunregister_fanin(SolverMatrix            * solvmtx,
                         starpu_data_handle_t  *** Lfanin_handle,
                         starpu_data_handle_t  *** Ufanin_handle);

int
starpu_zunregister_halo(SolverMatrix          * datacode,
                        starpu_data_handle_t ** Lhalo_handle,
                        starpu_data_handle_t ** Uhalo_handle);

int
starpu_zunregister_cblk( SolverMatrix          * datacode,
                         starpu_data_handle_t ** L_handle,
                         starpu_data_handle_t ** U_handle );

int
starpu_zunregister_data( Sopalin_Data_t * sopalin_data,
                         starpu_data_handle_t ** L_handle,
                         starpu_data_handle_t ** U_handle,
                         starpu_data_handle_t ** Lhalo_handle,
                         starpu_data_handle_t ** Uhalo_handle,
                         starpu_data_handle_t *** Lfanin_handle,
                         starpu_data_handle_t *** Ufanin_handle);
