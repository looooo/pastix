/**
 * @file drivers.h
 *
 *  $COPYRIGHTS$
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#ifndef _DRIVER_H_
#define _DRIVER_H_

#include "csc.h"

typedef enum pastix_driver_e {
    PastixDriverRSA, /* ok */
    PastixDriverCCC,//
    PastixDriverRCC,//
    PastixDriverOlaf,//
    PastixDriverPeer,//
    PastixDriverHB, /* ok */
    PastixDriverIJV, /* ok */
    PastixDriverMM, /* ok */
    PastixDriverDMM, /* ok */
    PastixDriverPetscS, /* ok */
    PastixDriverPetscU, /* ok */
    PastixDriverPetscH, /* ok */
    PastixDriverCSCD,//
    PastixDriverLaplacian, /* ok */
    PastixDriverBRGM,//
    PastixDriverBRGMD,//
    PastixDriverGraph
} pastix_driver_t;

void pastix_ex_getoptions(int argc, char **argv,
                          pastix_int_t *iparam, double *dparam,
                          pastix_driver_t *driver, char **filename );


void convertArrayToComplex64( pastix_int_t n, const double *A, void **B );
void convertArrayToComplex32( pastix_int_t n, const double *A, void **B );
void convertArrayToDouble(    pastix_int_t n, const double *A, void **B );
void convertArrayToFloat(     pastix_int_t n, const double *A, void **B );

int cscReadFromFile( pastix_driver_t  driver,
                     char            *filename,
                     pastix_csc_t    *csc,
                     MPI_Comm         pastix_comm );

int readHB   ( const char *filename, pastix_csc_t *csc );
int readRSA  ( const char *filename, pastix_csc_t *csc );
int readIJV  ( const char *filename, pastix_csc_t *csc );
int readMM   ( const char *filename, pastix_csc_t *csc );
int readDMM  ( const char *filename, pastix_csc_t *csc );
int readPETSC( const char *filename, pastix_csc_t *csc );
int readCSCD ( const char *filename, pastix_csc_t *csc, void **rhs, MPI_Comm pastix_comm );
int genLaplacian( const char *filename, pastix_csc_t *csc );

#endif /* _DRIVER_H_ */
