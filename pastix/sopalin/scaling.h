#ifndef SCALING_H
#define SCALING_H

void Matrix_Unscale_Sym(pastix_data_t * pastix_data, d_SolverMatrix *solvmtx, pastix_float_t *scaletab, pastix_float_t *iscaletab);

void Matrix_Unscale_Unsym(pastix_data_t * pastix_data, d_SolverMatrix *solvmtx, pastix_float_t *scalerowtab, pastix_float_t *iscalerowtab, pastix_float_t *scalecoltab, pastix_float_t *iscalecoltab);

#endif
