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
 *
 **/
#ifndef _DRIVER_H_
#define _DRIVER_H_

#include "csc.h"

typedef enum pastix_driver_e {
    PastixDriverRSA,
    PastixDriverCHB,
    PastixDriverCCC,
    PastixDriverRCC,
    PastixDriverOlaf,
    PastixDriverPeer,
    PastixDriverHB,
    PastixDriverIJV,
    PastixDriverMM,
    PastixDriverDMM,
    PastixDriverPetscS,
    PastixDriverPetscU,
    PastixDriverPetscH,
    PastixDriverCSCD,
    PastixDriverLaplacian,
    PastixDriverBRGM,
    PastixDriverBRGMD
} pastix_driver_t;

void pastix_ex_getoptions(int argc, char **argv,
                          pastix_int_t *iparam, double *dparam,
                          pastix_driver_t *driver, char **filename);


void convertArrayToComplex64( pastix_int_t n, const double *A, void **B );
void convertArrayToComplex32( pastix_int_t n, const double *A, void **B );
void convertArrayToDouble(    pastix_int_t n, const double *A, void **B );
void convertArrayToFloat(     pastix_int_t n, const double *A, void **B );

int cscReadFromFile( pastix_driver_t  driver,
                     char            *filename,
                     pastix_csc_t    *csc,
                     MPI_Comm         pastix_comm );


#endif /* _DRIVER_H_ */
