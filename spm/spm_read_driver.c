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
#include "spm.h"
#include "spm_drivers.h"
#if defined(HAVE_SCOTCH)
#include <scotch.h>
#endif

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Import a matrix file into a spm structure
 *
 * This function read or generate a sparse matrix from a file to store it into a
 * spm structure. The different formats accepted by this driver are described by
 * the driver field.
 *
 *******************************************************************************
 *
 * @param[in] driver
 *          This defines the driver to use to create the spm structure.
 *          = PastixDriverRSA
 *          = PastixDriverHB
 *          = PastixDriverIJV
 *          = PastixDriverMM
 *          = PastixDriverLaplacian
 *          = PastixDriverXLaplacian
 *
 * @param[in] filename
 *          The name of the file that stores the matrix (see driver)
 *
 * @param[in,out] spm
 *          On entry, an allocated sparse matrix structure.
 *          On exit, the filled sparse matrix structure with the matrix from the
 *          file.
 *
 * @param[in] pastix_comm
 *          The MPI communicator on which to distribute the sparse matrix. This
 *          is also used in case of distributed formats.
 *
 *******************************************************************************/
int
spmReadDriver( pastix_driver_t  driver,
               char            *filename,
               pastix_spm_t    *spm,
               MPI_Comm         comm )
{
    int mpirank = 0;
    int mpiinit;

    spm->mtxtype  = PastixGeneral;
    spm->flttype  = PastixDouble;
    spm->fmttype  = PastixCSC;
    spm->gN       = 0;
    spm->n        = 0;
    spm->gnnz     = 0;
    spm->nnz      = 0;
    spm->dof      = 1;
    spm->colptr   = NULL;
    spm->rowptr   = NULL;
    spm->values   = NULL;
    spm->loc2glob = NULL;

    MPI_Initialized( &mpiinit );
    if (mpiinit) {
        MPI_Comm_rank( comm, &mpirank );
    }

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
            readHB( filename, spm );
            break;

        case PastixDriverIJV:
            printf("driver: 3files file: %s\n", filename);
            readIJV( filename, spm );
            break;

        case PastixDriverMM:
            printf("driver: MatrixMarket file: %s\n", filename);
            readMM( filename, spm );
            break;

        case PastixDriverDMM:
            printf("driver: DistributedMatrixMarket file: %s\n", filename);
            //readMMD( filename, spm );
            break;

        case PastixDriverPetscS:
        case PastixDriverPetscU:
        case PastixDriverPetscH:
            printf("driver: PETSc file: %s\n", filename);
            //readPETSC( filename, spm );
            if (driver == PastixDriverPetscS) spm->mtxtype = PastixSymmetric;
            if (driver == PastixDriverPetscH) spm->mtxtype = PastixHermitian;
            break;

        case PastixDriverCSCD:
            printf("driver CSCd file: %s\n", filename);
            //readCSCD( filename, spm, rhs, comm );
            break;

        case PastixDriverLaplacian:
            if (mpirank == 0)
                printf("driver Laplacian: %s\n", filename);
            genLaplacian( filename, spm );
            break;

        case PastixDriverXLaplacian:
            if (mpirank == 0)
                printf("driver Extended Laplacian: %s\n", filename);
            genExtendedLaplacian( filename, spm );
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
            SCOTCH_graphData( &sgraph, NULL, &(spm->n), &(spm->colptr), NULL, NULL, NULL, NULL, &(spm->rowptr), NULL );
            fclose(file);

            spm->flttype = PastixPattern;
            spm->gN   = spm->n;
            spm->gnnz = spm->colptr[ spm->n ];
            spm->nnz  = spm->gnnz;
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
            readRSA( filename, spm );
        }

        spmConvert( PastixCSC, spm );
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
    /*       MPI_Comm_size(comm,&comm_size); */
    /*       send[0] = *ncol; */
    /*       send[1] = 0; */
    /*       if (*ncol != 0 && *rhs == NULL) */
    /*         send[1] = 1; */
    /*       MPI_Allreduce(send, recv, 2, PASTIX_MPI_INT, MPI_SUM, comm); */
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
    /*           MPI_Allreduce(RHS, RHS_recv, N, PASTIX_MPI_FLOAT, MPI_SUM, comm); */
    /*           free(RHS); */
    /*           *rhs = (pastix_complex64_t *) malloc((*ncol)*sizeof(pastix_complex64_t)); */

    /*           for (i = 0; i < *ncol; i++) */
    /*             (*rhs)[i] = RHS_recv[(*loc2glob)[i]-1]; */
    /*         } */
    /*     } */

    if ( mpiinit )
    {
        pastix_int_t nnz;

        if (mpirank == 0) {
            nnz = spm->colptr[spm->gN] - spm->colptr[0];
        }

        MPI_Bcast( spm, 2*sizeof(int)+3*sizeof(pastix_int_t), MPI_CHAR, 0, comm );
        MPI_Bcast( &nnz, 1, PASTIX_MPI_INT, 0, comm );

        fprintf(stderr, "%d: mtxtype=%d, flttype=%d, nnz=%ld, gN=%ld\n",
                mpirank, spm->mtxtype, spm->flttype, (long)nnz, (long)spm->gN );

        if ( mpirank != 0 )
        {
            spm->colptr = (pastix_int_t *) malloc((spm->gN+1) * sizeof(pastix_int_t));
            spm->rowptr = (pastix_int_t *) malloc(nnz * sizeof(pastix_int_t));
            spm->values = (void *)         malloc(nnz * pastix_size_of( spm->flttype ));
            spm->loc2glob = NULL;
            spm->loc2glob = NULL;
            /* spm->rhs    = (void *) malloc((*ncol)*sizeof(pastix_complex64_t)); */
            /* spm->type   = (char *) malloc(4*sizeof(char)); */
        }

        MPI_Bcast( spm->colptr, spm->gN+1, PASTIX_MPI_INT, 0, comm );
        MPI_Bcast( spm->rowptr, nnz,       PASTIX_MPI_INT, 0, comm );
        MPI_Bcast( spm->values, nnz * pastix_size_of( spm->flttype ), MPI_CHAR, 0, comm );
        /* MPI_Bcast(*rhs,    *ncol,   PASTIX_MPI_FLOAT, 0, comm); */
        /* MPI_Bcast(*type,    4,      MPI_CHAR,         0, comm); */
    }

    return PASTIX_SUCCESS;
}
