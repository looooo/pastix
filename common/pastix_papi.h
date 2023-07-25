/**
 *
 * @file pastix_papi.h
 *
 * @copyright 2004-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * PaStiX memory tracking function.
 *
 * @version 6.3.0
 * @author Mohamed Aymane Kherraz
 * @date 2023-08-01
 *
 */
#ifndef _pastix_papi_h_
#define _pastix_papi_h_

#if defined(PASTIX_WITH_PAPI)
#include <papi.h>

int    papiEnergyInit();
void   papiEnergyStart();
double papiEnergyStop();
void   papiEnergyFinalize();

#else

static inline int    papiEnergyInit() { return 0; }
static inline void   papiEnergyStart() {}
static inline double papiEnergyStop() { return 0.; }
static inline void   papiEnergyFinalize() {}

#endif

#endif /* _pastix_papi_h_ */
