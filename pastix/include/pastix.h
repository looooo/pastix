/*
 *  File: pastix.h
 *
 *  Part of a parallel direct block solver.
 *
 *  These lines are the common data
 *  declarations for all modules.
 *
 *  Authors:
 *    Mathieu  Faverge    - faverge@labri.fr
 *    David    GOUDIN     - .
 *    Pascal   HENON      - henon@labri.fr
 *    Xavier   LACOSTE    - lacoste@labri.fr
 *    Francois PELLEGRINI - .
 *    Pierre   RAMET      - ramet@labri.fr
 *
 *  Dates:
 *    Version 0.0 - from 08 may 1998
 *                  to   08 jan 2001
 *    Version 1.0 - from 06 jun 2002
 *                  to   06 jun 2002
 */
#ifndef _PASTIX_H_
#define _PASTIX_H_

#include "pastix/config.h"
#include "pastix/api.h"
#include <math.h>
#if defined(HAVE_MPI)
#include <mpi.h>
#else
#include "pastix/nompi.h"
#endif
#include "pastix/datatypes.h"
#include "pastix/retro.h"

struct pastix_data_s;
typedef struct pastix_data_s pastix_data_t;

/** ****************************************************************************
 *
 *  PaStiX constants - Compatible with CBLAS & LAPACK
 *  The naming and numbering is consistent with:
 *
 *    1) CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz),
 *    2) C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/).
 *
 **/
#define PastixNoTrans       111
#define PastixTrans         112
#define PastixConjTrans     113
#define PastixGeneral       111
#define PastixSymmetric     112
#define PastixHermitian     113

#endif /* _PASTIX_H_ */
