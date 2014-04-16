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
int starpu_zsubmit_outgoing_fanin(Sopalin_Data_t * sopalin_data,
                                  SolverCblk     * fcblk,
                                  SolverCblk     * hcblk);
#define starpu_submit_outgoing_fanin starpu_zsubmit_outgoing_fanin
