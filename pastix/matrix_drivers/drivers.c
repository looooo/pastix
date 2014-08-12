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

#include "common.h"
#include "drivers.h"

void convertArrayToComplex64( pastix_int_t n,
                              const double *A,
                              void **B )
{
    pastix_int_t i;
    const double *lA = A;
    pastix_complex64_t *lB;

    (*B) = malloc (n * sizeof(pastix_complex64_t) );
    assert(*B);
    lB = (pastix_complex64_t*)(*B);

    for(i=0; i<n; i++, lA++, lB++) {
        *lB = (double)(*lA) + I*0.;
    }
}

void convertArrayToComplex32( pastix_int_t n,
                              const double *A,
                              void **B )
{
    pastix_int_t i;
    const double *lA = A;
    pastix_complex32_t *lB;

    (*B) = malloc (n * sizeof(pastix_complex32_t) );
    assert(*B);
    lB = (pastix_complex32_t*)(*B);

    for(i=0; i<n; i++, lA++, lB++) {
        *lB = (float)(*lA) + I*0.;
    }
}

void convertArrayToDouble( pastix_int_t n,
                           const double *A,
                           void **B )
{
    pastix_int_t i;

    *B = malloc ( n * sizeof(double) );
    assert(*B);

    memcpy( *B, A, n * sizeof(double) );
}

void convertArrayToFloat( pastix_int_t n,
                          const double *A,
                          void **B )
{
    pastix_int_t i;
    const double *lA = A;
    float *lB;

    *B = malloc (n * sizeof(float) );
    assert(*B);
    lB = (float*)(*B);

    for(i=0; i<n; i++, lA++, lB++) {
        *lB = (float)(*lA);
    }
}

int cscReadFromFile( pastix_driver_t  driver,
                     char            *filename,
                     pastix_csc_t    *csc,
                     MPI_Comm         pastix_comm )
{
    pastix_int_t nrows;
    pastix_int_t nnz;
    int mpirank;

    MPI_Comm_rank( pastix_comm, &mpirank );

    if ( (mpirank == 0)                    ||
         (driver == PastixDriverLaplacian) ||
         (driver == PastixDriverCSCD)      ||
         (driver == PastixDriverBRGMD)     ||
         (driver == PastixDriverDMM)         )
    {
        switch(driver)
        {
        case PastixDriverCHB:
          /* printf("driver: CHB file: %s\n", filename); */
          /* z_chbRead(filename, */
          /*         ncol, &nrows, &nnz, */
          /*         colptr, rows, values, */
          /*         type, rhstype, rhs); */
          break;
        case PastixDriverCCC:
          /* printf("driver: CCC file: %s\n", filename); */
          /* z_cccRead(filename, */
          /*         ncol, &nrows, &nnz, */
          /*         colptr, rows, values, */
          /*         type, rhstype); */
          break;
        case PastixDriverRCC:
          /* printf("driver: RCC file: %s\n", filename); */
          /* z_cccRead(filename, */
          /*         ncol, &nrows, &nnz, */
          /*         colptr, rows, values, */
          /*         type, rhstype); */
          /* (*type)[0]='R'; */
          break;
        case PastixDriverOlaf:
          /* printf("driver: OLAF file: %s\n", filename); */
          /* z_olafRead(filename, */
          /*          ncol, &nrows, &nnz, */
          /*          colptr, rows, values, */
          /*          type, rhstype, */
          /*          rhs); */
          break;
        case PastixDriverPeer:
          /* printf("driver: PEER file: %s\n", filename); */
          /* z_peerRead(filename, */
          /*          ncol, &nrows, &nnz, */
          /*          colptr, rows, values, */
          /*          type, rhstype, */
          /*          rhs); */
      /* peerRead2("rsaname", ncol, &nrows, &nnz, colptr, rows, values, type, rhstype, rhs); */
          break;
        case PastixDriverHB:
          /* printf("driver: HB file: %s\n", filename); */
          /* z_HBRead(filename, */
          /*        ncol, &nrows, &nnz, */
          /*        colptr, rows, values, */
          /*        type, rhstype); */
          break;
        case PastixDriverIJV:
          /* printf("driver: 3files file: %s\n", filename); */
          /* z_threeFilesRead(filename, */
          /*                ncol, &nrows, &nnz, */
          /*                colptr, rows, values, */
          /*                type, rhstype); */
          break;
        case PastixDriverMM:
          /* printf("driver: MatrixMarket file: %s\n", filename); */
          /* z_MatrixMarketRead(filename, */
          /*                  ncol, &nrows, &nnz, */
          /*                  colptr, rows, values, */
          /*                  type, rhstype); */
          break;
        case PastixDriverDMM:
          /* printf("driver: DistributedMatrixMarket file: %s\n", filename); */
          /* z_DistributedMatrixMarketRead(filename, */
          /*                             ncol, &nrows, &nnz, */
          /*                             colptr, rows, values, loc2glob, */
          /*                             type, rhstype); */
          break;
        case PastixDriverPetscS:
        case PastixDriverPetscU:
        case PastixDriverPetscH:
          /* printf("driver: PETSc file: %s\n", filename); */
          /* z_PETScRead(filename, */
          /*           ncol, &nrows, &nnz, */
          /*           colptr, rows, values, */
          /*           type, rhstype); */
          /* if (driver_type == PETSCS) *type[1] = 'S'; */
          /* if (driver_type == PETSCH) *type[1] = 'H'; */
          break;
        case PastixDriverCSCD:
          /* printf("driver CSCdt file: %s\n", filename); */
          /* z_cscdRead(filename, */
          /*          colptr, rows, loc2glob, values, */
          /*          rhs, ncol, &nnz, */
          /*          pastix_comm); */

          /* *type = (char *) malloc(4*sizeof(char)); */
          /* sprintf(*type,"RSA"); */
          break;
        case PastixDriverLaplacian:
          /* if (mpirank == 0) */
          /*   printf("driver Laplacian\n"); */
          /* z_genlaplacian(*ncol, &nnz, colptr, rows, values, rhs, type, rhstype); */
          /* return EXIT_SUCCESS; */
          break;
/* #ifdef FDUPROS */
/*         case FDUP: */
/*           printf("driver: FDupros file: %s\n",filename); */
/*           z_driverFdupros(filename, */
/*                         ncol, &nrows, &nnz, */
/*                         colptr, rows, values, */
/*                         rhs, */
/*                         type, rhstype); */
/*           break; */
/*         case FDUP_DIST: */
/*           printf("driver: FDupros file: %s\n",filename); */
/*           z_driverFdupros_dist(filename, */
/*                              ncol, &nrows, &nnz, */
/*                              colptr, rows, loc2glob, values, */
/*                              rhs, */
/*                              type, rhstype, */
/*                              pastix_comm); */
/*           free(*rhs); */
/*           *rhs = NULL; */
/*           break; */
/* #endif */
        case PastixDriverRSA:
        default:
            rsaRead( filename, csc );
        }
    }

/* #ifndef TYPE_COMPLEX */
/*   if (*type) */
/*     if ((*type)[1] == 'H') */
/*       (*type)[1] = 'S'; */
/* #endif */

/*     /\* read RHS file *\/ */
/*   if (mpirank == 0) */
/*     { */
/*       FILE *file; */
/*       char filename2[256]; */
/*       pastix_int_t i; */
/*       double re; */

/*       sprintf(filename2,"%s.rhs",filename); */
/*       fprintf(stderr,"open RHS file : %s\n",filename2); */
/*       file = fopen(filename2,"r"); */
/*       if (file==NULL) */
/*         { */
/*           fprintf(stderr,"cannot load %s\n", filename2); */
/*         } */
/*       else */
/*         { */
/*           *rhs = (pastix_complex64_t *) malloc((*ncol)*sizeof(pastix_complex64_t)); */
/*           for (i = 0; i < *ncol; i++) */
/*             { */
/*               (*rhs)[i] = 0.0; */
/*               if (1 != fscanf(file,"%lg\n", &re)) */
/*                 { */
/*                   fprintf(stderr, "ERROR: reading rhs(%ld)\n", (long int)i); */
/*                   exit(1); */
/*                 } */
/*               (*rhs)[i] = (pastix_complex64_t)re; */
/*             } */
/*           fclose(file); */
/*         } */
/*     } */

/*   if (*rhs == NULL  && ( mpirank == 0 || */
/*                          driver_type == LAPLACIAN) && */
/*       driver_type != CSCD && driver_type != FDUP_DIST && driver_type != MMD) */
/*     { */
/*       *rhs = (pastix_complex64_t *) malloc((*ncol)*sizeof(pastix_complex64_t)); */
/*       pastix_int_t i,j; */
/*       for (i = 0; i < *ncol; i++) */
/*         (*rhs)[i] = 0.0; */

/* #ifdef TYPE_COMPLEX */
/*       fprintf(stdout, "Setting right-hand-side member such as X[i] = i + i*I\n"); */
/* #else */
/*       fprintf(stdout, "Setting right-hand-side member such as X[i] = i\n"); */
/* #endif */
/*       for (i = 0; i < *ncol; i++) */
/*         { */
/*           for (j = (*colptr)[i] -1; j < (*colptr)[i+1]-1; j++) */
/*             { */
/* #ifdef TYPE_COMPLEX */
/*               (*rhs)[(*rows)[j]-1] += (i+1 + I*(i+1))*(*values)[j]; */
/* #else */
/*               (*rhs)[(*rows)[j]-1] += (i+1)*(*values)[j]; */
/* #endif */
/*               if (MTX_ISSYM((*type)) && i != (*rows)[j]-1) */
/*                 { */
/* #ifdef TYPE_COMPLEX */
/*                   (*rhs)[i] += ((*rows)[j] + (*rows)[j] * I)*(*values)[j]; */
/* #else */
/*                   (*rhs)[i] += ((*rows)[j])*(*values)[j]; */
/* #endif */
/*                 } */
/*             } */
/*         } */
/*     } */
/*   if (driver_type == CSCD || driver_type == FDUP_DIST || driver_type == MMD) */
/*     { */
/*       pastix_int_t i,j; */
/*       pastix_int_t N; */
/*       pastix_int_t send[2], recv[2]; */
/*       int comm_size; */
/*       MPI_Comm_size(pastix_comm,&comm_size); */
/*       send[0] = *ncol; */
/*       send[1] = 0; */
/*       if (*ncol != 0 && *rhs == NULL) */
/*         send[1] = 1; */
/*       MPI_Allreduce(send, recv, 2, PASTIX_MPI_INT, MPI_SUM, pastix_comm); */
/*       N = recv[0]; */
/*       if (recv[1] > 0) */
/*         { */
/*           pastix_complex64_t *RHS; */
/*           pastix_complex64_t *RHS_recv; */
/*           RHS  = (pastix_complex64_t *) malloc((N)*sizeof(pastix_complex64_t)); */

/*           for (i = 0; i < N; i++) */
/*             (RHS)[i] = 0.0; */
/*           if (mpirank == 0) */
/*             fprintf(stdout, "Setting RHS such as X[i] = i\n"); */
/*           for (i = 0; i < *ncol; i++) */
/*             { */
/*               for (j = (*colptr)[i] -1; j < (*colptr)[i+1]-1; j++) */
/*                 { */
/*                   (RHS)[(*rows)[j]-1] += ((*loc2glob)[i])*(*values)[j]; */
/*                   if (MTX_ISSYM((*type)) && i != (*rows)[j]-1) */
/*                     (RHS)[(*loc2glob)[i]-1] += ((*rows)[j])*(*values)[j]; */
/*                 } */
/*             } */
/*           RHS_recv  = (pastix_complex64_t *) malloc((N)*sizeof(pastix_complex64_t)); */
/*           MPI_Allreduce(RHS, RHS_recv, N, PASTIX_MPI_FLOAT, MPI_SUM, pastix_comm); */
/*           free(RHS); */
/*           *rhs = (pastix_complex64_t *) malloc((*ncol)*sizeof(pastix_complex64_t)); */

/*           for (i = 0; i < *ncol; i++) */
/*             (*rhs)[i] = RHS_recv[(*loc2glob)[i]-1]; */
/*         } */
/*     } */

/*   if ( driver_type != LAPLACIAN && driver_type != CSCD && */
/*        driver_type != FDUP_DIST && driver_type != MMD) */
/*     { */
/*       MPI_Bcast(ncol,1,PASTIX_MPI_INT,0,pastix_comm); */
/*       MPI_Bcast(&nnz,1,PASTIX_MPI_INT,0,pastix_comm); */

/*       if (mpirank!=0) */
/*         { */
/*           *colptr = (pastix_int_t *)   malloc((*ncol+1)*sizeof(pastix_int_t)); */
/*           *rows   = (pastix_int_t *)   malloc(nnz*sizeof(pastix_int_t)); */
/*           *values = (pastix_complex64_t *) malloc(nnz*sizeof(pastix_complex64_t)); */
/*           *rhs    = (pastix_complex64_t *) malloc((*ncol)*sizeof(pastix_complex64_t)); */
/*           *type   = (char *)  malloc(4*sizeof(char)); */
/*         } */

/*       MPI_Bcast(*colptr, *ncol+1, PASTIX_MPI_INT,   0, pastix_comm); */
/*       MPI_Bcast(*rows,    nnz,    PASTIX_MPI_INT,   0, pastix_comm); */
/*       MPI_Bcast(*values,  nnz,    PASTIX_MPI_FLOAT, 0, pastix_comm); */
/*       MPI_Bcast(*rhs,    *ncol,   PASTIX_MPI_FLOAT, 0, pastix_comm); */
/*       MPI_Bcast(*type,    4,      MPI_CHAR,         0, pastix_comm); */
/*   } */

  return EXIT_SUCCESS;
}

