#ifndef SCALING_H
#define SCALING_H

void Matrix_Unscale_Sym(pastix_data_t * pastix_data, SolverMatrix *solvmtx, PASTIX_FLOAT *scaletab, PASTIX_FLOAT *iscaletab);

void Matrix_Unscale_Unsym(pastix_data_t * pastix_data, SolverMatrix *solvmtx, PASTIX_FLOAT *scalerowtab, PASTIX_FLOAT *iscalerowtab, PASTIX_FLOAT *scalecoltab, PASTIX_FLOAT *iscalecoltab);

#endif
