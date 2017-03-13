/**
 * @file spm_drivers.h
 *
 * SParse Matrix package driver header.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#ifndef _SPM_DRIVER_H_
#define _SPM_DRIVER_H_

#include "spm.h"

void convertArrayToComplex64( pastix_int_t n, const double *A, void **B );
void convertArrayToComplex32( pastix_int_t n, const double *A, void **B );
void convertArrayToDouble(    pastix_int_t n, const double *A, void **B );
void convertArrayToFloat(     pastix_int_t n, const double *A, void **B );

int readHB   ( const char *filename, pastix_spm_t *spm );
int readRSA  ( const char *filename, pastix_spm_t *spm );
int readIJV  ( const char *filename, pastix_spm_t *spm );
int readMM   ( const char *filename, pastix_spm_t *spm );
int readDMM  ( const char *filename, pastix_spm_t *spm );
int readPETSC( const char *filename, pastix_spm_t *spm );
int readCSCD ( const char *filename, pastix_spm_t *spm, void **rhs, MPI_Comm pastix_comm );
int genLaplacian( const char *filename, pastix_spm_t *spm );
int genExtendedLaplacian( const char *filename, pastix_spm_t *spm );

#endif /* _SPM_DRIVER_H_ */
