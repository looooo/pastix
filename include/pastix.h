/**
 *
 * @file pastix.h
 *
 * PaStiX main header file.
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Alban Bellot
 * @author Astrid Casadei
 * @author David Goudin
 * @author Francois Pellegrini
 * @author Gregoire Pichon
 * @author Mathias Hastaran
 * @author Mathieu Faverge
 * @author Pascal Henon
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#ifndef _PASTIX_H_
#define _PASTIX_H_

#include "pastix/config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include "pastix/api.h"
#include "pastix/datatypes.h"

#include <math.h>
#if defined(HAVE_MPI)
#include <mpi.h>
#else
#include "pastix/nompi.h"
#endif

/*
 * Main function for compatibility with former versions
 */
int pastix( pastix_data_t **pastix_data,
            MPI_Comm        pastix_comm,
            pastix_int_t    n,
            pastix_int_t   *colptr,
            pastix_int_t   *row,
            void           *avals,
            pastix_int_t   *perm,
            pastix_int_t   *invp,
            void           *b,
            pastix_int_t    nrhs,
            pastix_int_t   *iparm,
            double         *dparm );

/*
 * Solver initialization
 */
void pastixInitParam( pastix_int_t   *iparm,
                      double         *dparm );
void pastixInit     ( pastix_data_t **pastix_data,
                      MPI_Comm        pastix_comm,
                      pastix_int_t   *iparm,
                      double         *dparm );
void pastixFinalize ( pastix_data_t **pastix_data );

/*
 * Main steps of the solver
 */
int pastix_task_analyze( pastix_data_t      *pastix_data,
                         pastix_spm_t       *spm );
int pastix_task_numfact( pastix_data_t      *pastix_data,
                         pastix_spm_t       *spm );
int pastix_task_solve  ( pastix_data_t      *pastix_data,
                         const pastix_spm_t *spm,
                         pastix_int_t        nrhs,
                         void               *b,
                         pastix_int_t        ldb );
int pastix_task_refine ( pastix_data_t      *pastix_data,
                         void               *x,
                         pastix_int_t        rhsnbr,
                         void               *b );

/*
 * Analyze subtasks
 */
int pastix_subtask_order     ( pastix_data_t      *pastix_data,
                               const pastix_spm_t *spm,
                               Order              *myorder );
int pastix_subtask_symbfact  ( pastix_data_t      *pastix_data );
int pastix_subtask_reordering( pastix_data_t      *pastix_data );
int pastix_subtask_blend     ( pastix_data_t      *pastix_data );

/*
 * Numerical factorization subtasks
 */
int pastix_subtask_spm2bcsc  ( pastix_data_t      *pastix_data,
                               pastix_spm_t       *spm );
int pastix_subtask_bcsc2ctab ( pastix_data_t      *pastix_data,
                               const pastix_spm_t *spm );
int pastix_subtask_sopalin   ( pastix_data_t      *pastix_data,
                               const pastix_spm_t *spm );

/*
 * Numerical solve subtasks
 */
int pastix_subtask_applyorder( pastix_data_t    *pastix_data,
                               pastix_coeftype_t flttype,
                               pastix_dir_t      dir,
                               pastix_int_t      m,
                               pastix_int_t      n,
                               void             *b,
                               pastix_int_t      ldb );
int pastix_subtask_trsm( pastix_data_t    *pastix_data,
                         pastix_coeftype_t flttype,
                         pastix_side_t     side,
                         pastix_uplo_t     uplo,
                         pastix_trans_t    trans,
                         pastix_diag_t     diag,
                         pastix_int_t      nrhs,
                         void             *b,
                         pastix_int_t      ldb );
int pastix_subtask_diag( pastix_data_t    *pastix_data,
                         pastix_coeftype_t flttype,
                         pastix_int_t      nrhs,
                         void             *b,
                         pastix_int_t      ldb );

/*
 * Schur complement manipulation routines.
 */
void pastix_setSchurUnknownList( pastix_data_t       *pastix_data,
                                 pastix_int_t         n,
                                 const pastix_int_t  *list );
int  pastix_getSchur           ( const pastix_data_t *pastix_data,
                                 void                *S,
                                 pastix_int_t         lds );

#endif /* _PASTIX_H_ */
