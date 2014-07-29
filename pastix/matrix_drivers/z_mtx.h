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
#ifndef _MTX_H_
#define _MTX_H_

#include <stdio.h>

#ifndef MTX_ISSYM
#define MTX_ISSYM(a) ( (a)[1] == 'S'  )
#define MTX_ISHER(a) ( (a)[1] == 'H'  )
#define MTX_ISCOM(a) ( (a)[0] == 'C'  )
#define MTX_ISRHX(a) ( (a)[2] == 'X'  )
#define MTX_ISRHS(a) ( (a)[0] != '\0' )
#endif

void z_MatrixMarketRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_complex64_t **val, char **Type, char **RhsType);
void z_mumpsReadHeader(FILE *infile, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type);
void  z_mumpsRead(char const *dirname, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_complex64_t **val, char **Type, char **RhsType);
void mtxSym(pastix_int_t Nrow, pastix_int_t Ncol, pastix_int_t Nnzero, pastix_int_t *col, pastix_int_t *row, pastix_complex64_t *val);
void mtxCheck(pastix_int_t Nrow, pastix_int_t Ncol, pastix_int_t Nnzero, pastix_int_t *col, pastix_int_t *row, pastix_complex64_t *val);
void z_mtxReadHeader(FILE *infile, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type);
void z_mtxRead(char *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_complex64_t **val, char **Type, char **RhsType, pastix_int_t flagsort);
void z_ijvReadHeader(FILE *infile, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type);
void z_ijvRead(char *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_complex64_t **val, char **Type, char **RhsType, pastix_int_t flagsort);
void z_rsaReadHeader(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type, char *RhsType);
void z_rsaRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_complex64_t **val, char **Type, char **RhsType);
#ifdef PREC_DOUBLE
void z_HBRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_complex64_t **val, char **Type, char **RhsType);
#endif
void z_chbRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_complex64_t **val, char **Type, char **RhsType, pastix_complex64_t **rhs);
void z_cccReadHeader(FILE *infile, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type);
void z_cccRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_complex64_t **val, char **Type, char **RhsType);
void z_olafReadHeader(FILE *infile, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type);
void z_olafRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_complex64_t **val, char **Type, char **RhsType, pastix_complex64_t **rhs);
void z_peerRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col,pastix_int_t **row, pastix_complex64_t **val, char **Type, char **RhsType, pastix_complex64_t **rhs);

void z_diag_dominance(pastix_int_t baseval, pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_complex64_t *a);
void z_diag_unite(pastix_int_t baseval, pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_complex64_t *a);
void no_diag(pastix_int_t baseval, pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja);
void symmetrize_pattern(pastix_int_t n, pastix_int_t **iaptr, pastix_int_t **japtr, pastix_complex64_t **aaptr);
void z_dimsym(pastix_int_t n, pastix_int_t **iaptr, pastix_int_t **japtr);
void z_checkStrucSym(pastix_int_t n, pastix_int_t *nz, pastix_int_t **colptr, pastix_int_t **row, pastix_complex64_t **avals);

void z_driverFdupros(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_complex64_t **val,
		   pastix_complex64_t ** rhs, char **Type, char **RhsType);

#include "z_csparse.h"

#endif
