/**
 * @file pastix_papi.c
 *
 * @brief A file that contains all the functions used for papi manipulation.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Mohamed Aymane Kherraz
 * @date 2023-08-01
 *
 * @code
 *
 */
#include "common.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

static const char *_pastix_papi_event_names[] =
{
    "rapl::RAPL_ENERGY_PKG:cpu=%d",
    "rapl::RAPL_ENERGY_DRAM:cpu=%d"
};

static int _pastix_papi_N_EVTS      = 2;
static int _pastix_papi_N_SOCK      = 0;
static int _pastix_papi_initialized = 0;
static int _pastix_papi_EventSet    = PAPI_NULL;

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief A helper function used internally to add PAPI events to the event set.
 *
 *******************************************************************************
 *
 * @param[in] EventSet
 *          Represents the total time taken to complete a specific task or
 *          operation.
 *
 * @param[in] socket
 *          A pointer to a structure that stores energy consumption and average
 *          power consumption values.
 *
 *******************************************************************************
 *
 * @return PASTIX_SUCCESS on success, PASTIX_ERR_INTERNAL if any of the event
 *         could not be added to the EventSet.
 *
 *******************************************************************************/
static inline int
pastix_papi_add_event( int EventSet,
                       int socket )
{
    int i, retval, rc = PASTIX_SUCCESS;
    for ( i = 0; i < _pastix_papi_N_EVTS; i++ ) {
        char buf[256];
        int  code;

        PAPI_event_info_t info;

        snprintf( buf, 255, _pastix_papi_event_names[i], socket );
        buf[255] = '\0';
        retval = PAPI_event_name_to_code( buf, &code );
        if ( retval != PAPI_OK ) {
            pastix_print_warning( "Failed to convert PAPI %s event into code\n", buf );
            rc = PASTIX_ERR_INTERNAL;
            continue;
        }

        retval = PAPI_get_event_info( code, &info );
        if ( retval != PAPI_OK ) {
            pastix_print_warning( "Failed to get PAPI event information for event %s\n", buf );
            rc = PASTIX_ERR_INTERNAL;
            continue;
        }

        retval = PAPI_add_event( EventSet, code );
        if ( retval != PAPI_OK ) {
            pastix_print_warning( "Failed to add PAPI event %s to the eventset\n", buf );
            rc = PASTIX_ERR_INTERNAL;
        }
    }
    return rc;
}

/**
 *******************************************************************************
 *
 * @brief Init the PAPI energy counters.
 *
 *******************************************************************************
 *
 * @param[in] nbr_socks
 *          The number of sockets.
 *
 *******************************************************************************/
int
papiEnergyInit( pastix_int_t nbr_socks )
{
    int i, retval;

    /* Initialize the PAPI library */
    retval = PAPI_library_init( PAPI_VER_CURRENT );
    if ( retval != PAPI_VER_CURRENT ) {
        _pastix_papi_initialized = -1;
        pastix_print_warning( "Could not initialize PAPI\n");
        return PASTIX_ERR_INTERNAL;
    }

    /* Creating the eventset */
    retval = PAPI_create_eventset( &_pastix_papi_EventSet );
    if ( retval != PAPI_OK ) {
        _pastix_papi_initialized = -1;
        pastix_print_warning( "Failed to create the PAPI EventSet\n");
        return PASTIX_ERR_INTERNAL;
    }

    /* TODO: Check how to detect the number of sockets */
    _pastix_papi_N_SOCK = nbr_socks;

    /* Add the event code for each socket */
    for ( i = 0; i < _pastix_papi_N_SOCK; i++ ) {
        pastix_papi_add_event( _pastix_papi_EventSet, i );
    }

    _pastix_papi_initialized = 1;
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Starts counting the energy consumption from this point onwards.
 *
 *******************************************************************************/
void
papiEnergyStart()
{
    int retval;

    if ( _pastix_papi_initialized < 1 ) {
        return;
    }

    /* Start counting */
    retval = PAPI_start( _pastix_papi_EventSet );
    if ( retval != PAPI_OK ) {
        _pastix_papi_initialized = -1;
        pastix_print_warning( "Failed to start PAPI counters\n");
        return;
    }

    _pastix_papi_initialized = 2;
}

/**
 *******************************************************************************
 *
 * @brief Stops counting energy consumption and stores the energy consumption
 * values into an array.
 *
 *******************************************************************************
 *
 * @return The energy consumed on all sockets (sockets + DRAM) since the last
 * call to papiEnergyStart().
 *
 *******************************************************************************/
double
papiEnergyStop()
{
    long long energy = 0;
    long long values[_pastix_papi_N_EVTS * _pastix_papi_N_SOCK];
    int       retval, k;

    if ( _pastix_papi_initialized < 2 ) {
        if ( _pastix_papi_initialized > 0 ) {
            pastix_print_warning( "PAPI counters were not started\n");
        }
        return 0.0;
    }

    /* Stop counting and store the values into the array */
    retval = PAPI_stop( _pastix_papi_EventSet, values );
    if ( retval != PAPI_OK ) {
        pastix_print_warning( "Failed to stop PAPI counters\n");
        return 0.0;
    }

    _pastix_papi_initialized = 1;

    energy = 0;
    for ( k = 0; k < _pastix_papi_N_EVTS * _pastix_papi_N_SOCK; k++ ) {
        energy = energy + values[k];
    }

    return (double)energy;
}

/**
 *******************************************************************************
 *
 * @brief Finalizes the PAPI library and frees the resources used by PAPI after
 * the energy trace is complete.
 *
 *******************************************************************************/
void
papiEnergyFinalize( )
{
    /* free the resources used by PAPI */
    if ( _pastix_papi_initialized > 0 ) {
        PAPI_shutdown();
    }
    _pastix_papi_initialized = 0;
}
