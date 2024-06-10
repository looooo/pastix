/**
 *
 * @file pastix_papi.h
 *
 * @copyright 2004-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * PaStiX memory tracking function.
 *
 * @version 6.3.2
 * @author Mohamed Aymane Kherraz
 * @author Alycia Lisito
 * @author Mathieu Faverge
 * @date 2023-11-22
 *
 */
#ifndef _pastix_papi_h_
#define _pastix_papi_h_

#if defined(PASTIX_WITH_PAPI)

int    papiEnergyInit( pastix_int_t nbr_socks );
void   papiEnergyStart( void );
double papiEnergyStop( void );
void   papiEnergyFinalize( void );

#else

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static inline int    papiEnergyInit( pastix_int_t nbr_socks )
{
    (void)nbr_socks;
    return 0;
}
static inline void   papiEnergyStart( void ) {}
static inline double papiEnergyStop( void ) { return 0.; }
static inline void   papiEnergyFinalize( void ) {}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#endif

#endif /* _pastix_papi_h_ */
