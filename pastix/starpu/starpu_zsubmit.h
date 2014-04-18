/**
 * @file starpu_zsubmit.h
 *
 * @author Xavier Lacoste
 * @precisions normal z -> s d c
 */

#include "common.h"
#include "sopalin3d.h"
#include "solver.h"

int
starpu_zgesubmit_incomming_fanin(Sopalin_Data_t * sopalin_data);
int
starpu_zsysubmit_incomming_fanin(Sopalin_Data_t * sopalin_data);
int starpu_zgesubmit_outgoing_fanin(Sopalin_Data_t * sopalin_data,
                                    SolverCblk     * fcblk,
                                    SolverCblk     * hcblk);
int starpu_zsysubmit_outgoing_fanin(Sopalin_Data_t * sopalin_data,
                                    SolverCblk     * fcblk,
                                    SolverCblk     * hcblk);
#if defined(CHOL_SOPALIN) && defined(SOPALIN_LU)
#  define starpu_submit_outgoing_fanin starpu_zgesubmit_outgoing_fanin
#else
#  define starpu_submit_outgoing_fanin starpu_zsysubmit_outgoing_fanin
#endif
