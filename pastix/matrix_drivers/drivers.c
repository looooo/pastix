/**
 * @file drivers.c
 *
 *  $COPYRIGHTS$
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2011-11-11
 *
 **/

#include "common.h"
#include "drivers.h"
#if defined(HAVE_SCOTCH)
#include <scotch.h>
#endif

int cscReadFromFile( pastix_driver_t  driver,
                     char            *filename,
                     pastix_csc_t    *csc,
                     void           **rhs,
                     MPI_Comm         pastix_comm )
{
    int mpirank;

    csc->mtxtype  = PastixGeneral;
    csc->flttype  = PastixDouble;
    csc->fmttype  = PastixCSC;
    csc->gN       = 0;
    csc->n        = 0;
    csc->gnnz     = 0;
    csc->nnz      = 0;
    csc->dof      = 1;
    csc->colptr   = NULL;
    csc->rows     = NULL;
    csc->avals    = NULL;
    csc->loc2glob = NULL;

    MPI_Comm_rank( pastix_comm, &mpirank );

    if ( (mpirank == 0)                    ||
         (driver == PastixDriverLaplacian) ||
         (driver == PastixDriverCSCD)      ||
         (driver == PastixDriverBRGMD)     ||
         (driver == PastixDriverDMM)         )
    {
        switch(driver)
        {
        case PastixDriverCCC:
        case PastixDriverRCC:
        case PastixDriverOlaf:
        case PastixDriverPeer:
            fprintf(stderr, "driver: Driver not implemented\n");
            break;

        case PastixDriverHB:
            /* TODO: Possible to read the RHS, the solution or a guess of the solution */
            printf("driver: HB file: %s\n", filename);
            readHB( filename, csc );
            break;

        case PastixDriverIJV:
            printf("driver: 3files file: %s\n", filename);
            readIJV( filename, csc );
            break;

        case PastixDriverMM:
            printf("driver: MatrixMarket file: %s\n", filename);
            readMM( filename, csc );
            break;

        case PastixDriverDMM:
            printf("driver: DistributedMatrixMarket file: %s\n", filename);
            //readMMD( filename, csc );
            break;

        case PastixDriverPetscS:
        case PastixDriverPetscU:
        case PastixDriverPetscH:
            printf("driver: PETSc file: %s\n", filename);
            //readPETSC( filename, csc );
            if (driver == PastixDriverPetscS) csc->mtxtype = PastixSymmetric;
            if (driver == PastixDriverPetscH) csc->mtxtype = PastixHermitian;
            break;

        case PastixDriverCSCD:
            printf("driver CSCd file: %s\n", filename);
            //readCSCD( filename, csc, rhs, pastix_comm );
            break;

        case PastixDriverLaplacian:
            if (mpirank == 0)
                printf("driver Laplacian: %s\n", filename);
            genLaplacian( filename, csc );
            break;

        case PastixDriverGraph:
#if defined(HAVE_SCOTCH)
        {
            SCOTCH_Graph sgraph;
            FILE *file = fopen( filename, "r" );

            /* Check integer compatibility */
            if (sizeof(pastix_int_t) != sizeof(SCOTCH_Num)) {
                errorPrint("Inconsistent integer type\n");
                fclose(file);
                return PASTIX_ERR_INTEGER_TYPE;
            }

            SCOTCH_graphLoad( &sgraph, file, 1, 0 );
            SCOTCH_graphData( &sgraph, NULL, &(csc->n), &(csc->colptr), NULL, NULL, NULL, NULL, &(csc->rows), NULL );
            fclose(file);

            csc->flttype = PastixPattern;
            csc->gN   = csc->n;
            csc->gnnz = csc->colptr[ csc->n ];
            csc->nnz  = csc->gnnz;
        }
#else
        {
            fprintf(stderr, "Scotch driver to read graph file unavailable.\n"
                    "Compile with Scotch support to provide it\n");
            return PASTIX_ERR_BADPARAMETER;
        }
#endif
        break;

        case PastixDriverRSA:
        default:
            readRSA( filename, csc );
        }

        spmConvert( PastixCSC, csc );
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

    if ( 1 )
    {
        pastix_int_t nnz;

        fprintf(stderr, "Hello\n");
        if (mpirank == 0) {
            nnz = csc->colptr[csc->gN] - csc->colptr[0];
        }

        MPI_Bcast( csc, 2*sizeof(int)+3*sizeof(pastix_int_t), MPI_CHAR, 0, pastix_comm );
        MPI_Bcast( &nnz, 1, PASTIX_MPI_INT, 0, pastix_comm );

        fprintf(stderr, "%d: mtxtype=%d, flttype=%d, nnz=%ld, gN=%ld\n",
                mpirank, csc->mtxtype, csc->flttype, (long)nnz, (long)csc->gN );

        if ( mpirank != 0 )
        {
            csc->colptr = (pastix_int_t *) malloc((csc->gN+1) * sizeof(pastix_int_t));
            csc->rows   = (pastix_int_t *) malloc(nnz * sizeof(pastix_int_t));
            csc->avals  = (void *)         malloc(nnz * pastix_size_of( csc->flttype ));
            csc->loc2glob = NULL;
            csc->loc2glob = NULL;
            /* csc->rhs    = (void *) malloc((*ncol)*sizeof(pastix_complex64_t)); */
            /* csc->type   = (char *) malloc(4*sizeof(char)); */
        }

        MPI_Bcast( csc->colptr, csc->gN+1, PASTIX_MPI_INT, 0, pastix_comm );
        MPI_Bcast( csc->rows,   nnz,       PASTIX_MPI_INT, 0, pastix_comm );
        MPI_Bcast( csc->avals,  nnz * pastix_size_of( csc->flttype ), MPI_CHAR, 0, pastix_comm );
        /* MPI_Bcast(*rhs,    *ncol,   PASTIX_MPI_FLOAT, 0, pastix_comm); */
        /* MPI_Bcast(*type,    4,      MPI_CHAR,         0, pastix_comm); */
    }

    return PASTIX_SUCCESS;
}
