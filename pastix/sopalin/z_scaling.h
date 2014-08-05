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
#ifndef SCALING_H
#define SCALING_H

void z_Matrix_Unscale_Sym(z_pastix_data_t * pastix_data, z_SolverMatrix *solvmtx, pastix_complex64_t *scaletab, pastix_complex64_t *iscaletab);

void z_Matrix_Unscale_Unsym(z_pastix_data_t * pastix_data, z_SolverMatrix *solvmtx, pastix_complex64_t *scalerowtab, pastix_complex64_t *iscalerowtab, pastix_complex64_t *scalecoltab, pastix_complex64_t *iscalecoltab);

#endif
