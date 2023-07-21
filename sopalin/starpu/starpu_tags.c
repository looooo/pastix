/**
 *
 * @file starpu_tags.c
 *
 * @copyright 2017-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2023-01-17
 *
 * Functions to manage the MPI data tags.
 *
 * @addtogroup pastix_starpu
 * @{
 *
 **/
#include "common.h"
#if !defined(PASTIX_WITH_STARPU)
#error "This file should not be compiled if Starpu is not enabled"
#endif
#include "pastix_starpu.h"

#if defined(PASTIX_WITH_MPI)

/**
 * @brief Structure Pastix StarPU tag
 *
 * List structure to manage the set of available tags.
 */
struct pst_range_;
typedef struct pst_range_ pst_range_t;

struct pst_range_ {
    int64_t       min;  /**< Minimal value in the range     */
    int64_t       max;  /**< Maximal value in the range     */
    pst_range_t  *next; /**< Pointer to the following range */
};

/**
 * @brief Pointer to the first set or registered tags
 */
static pst_range_t *pst_first     = NULL;

/**
 * @brief StarPU tag upper bound
 */
static int64_t starpu_tag_ub = 0;

/**
 *******************************************************************************
 *
 * @brief Initialize the StarPU tags manager.
 *
 *******************************************************************************
 *
 * @param[in] pastix
 *          The main pastix_data structure to provide the MPI communicator.
 *
 *******************************************************************************
 *
 * @retval TODO
 *
 ******************************************************************************/
int
pastix_starpu_tag_init( pastix_data_t *pastix )
{
    if (!starpu_tag_ub) {
        int          ok       = 0;
        void        *tag_ub_p = NULL;

        starpu_mpi_comm_get_attr( pastix->inter_node_comm, STARPU_MPI_TAG_UB, &tag_ub_p, &ok );
        starpu_tag_ub = (uint64_t)((intptr_t)tag_ub_p);

        if ( !ok ) {
            pastix_print_error("pastix_starpu_tag_init: MPI_TAG_UB not known by StarPU\n");
        }

        return PASTIX_SUCCESS;
    }
    else {
        return PASTIX_ERR_INTERNAL;
    }
}

/**
 *******************************************************************************
 *
 * @brief Book a range of StarPU unique tags of size nbtags.
 *
 * This function returns the minimal tag value available to allow the
 * registration of nbtags data in a continuous range.
 *
 * Note that this function must be called exactly the same way on all nodes to
 * make sure the tags are identical from one node to another.
 *
 *******************************************************************************
 *
 * @param[in] nbtags
 *          The number of tags required to register the sparse matrix or right
 *          hand side.
 *
 *******************************************************************************
 *
 * @return V, the minimal tag value to use. The range [V:V+nbtags-1] is booked.
 *
 ********************************************************************************/
int64_t
pastix_starpu_tag_book( int64_t nbtags )
{
    pst_range_t *new;
    pst_range_t *prev    = NULL;
    pst_range_t *current = pst_first;
    int64_t      min = 0;
    int64_t      max = ( current == NULL ) ? starpu_tag_ub : current->min;

    assert( starpu_tag_ub != 0 ); /* StarPU tag must be initialized */

    while ( ((max - min) < nbtags) && (current != NULL) ) {
        min     = current->max;
        prev    = current;
        current = current->next;
        max     = ( current == NULL ) ? starpu_tag_ub : current->min;
    }

    if ( (max - min) < nbtags ) {
        pastix_print_error( "pastix_starpu_tag_book: No space left in tags (looking for %ld tags)\n",
                            nbtags );
        return -1;
    }

    new = malloc( sizeof( pst_range_t ) );
    new->min  = min;
    new->max  = min + nbtags;
    new->next = current;
    if ( prev == NULL ) {
        pst_first = new;
    }
    else {
        assert( prev->next == current );
        prev->next = new;
    }

#if defined(PASTIX_DEBUG_STARPU)
    fprintf( stderr, "pastix_starpu_tag: Book %ld - %ld\n",
             min, min + nbtags );
#endif

    assert( pst_first != NULL );
    return new->min;
}

/**
 *******************************************************************************
 *
 * @brief Release the set of tags starting by min.
 *
 * This function releases the range of tags that starts by the min value.
 *
 *******************************************************************************
 *
 * @param[in] min
 *          The initial value in the range
 *
 ******************************************************************************/
void
pastix_starpu_tag_release( int64_t min )
{
    pst_range_t *prev    = NULL;
    pst_range_t *current = pst_first;

    assert( pst_first != NULL ); /* At least one range must be registered */

    while ( (current != NULL) && (current->min < min) ) {
        prev    = current;
        current = current->next;
    }

    assert( current != NULL );
    assert( current->min == min );

    if ( prev ) {
        prev->next = current->next;
    }
    else {
        assert( current == pst_first );
        pst_first = current->next;
    }

#if defined(PASTIX_DEBUG_STARPU)
    fprintf( stderr, "pastix_starpu_tag: Release %ld - %ld\n",
             current->min, current->max );
#endif

    free( current );

    return;
}

#else /* defined(PASTIX_WITH_MPI) */

/**
 *******************************************************************************
 *
 * @brief Initialize the StarPU tags manager.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The main pastix_data structure to provide the MPI communicator.
 *
 *******************************************************************************
 *
 * @retval TODO
 *
 ******************************************************************************/
int
pastix_starpu_tag_init( __attribute__((unused)) pastix_data_t *pastix_data ) {
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Book a range of StarPU unique tags of size nbtags.
 *
 * This function returns the minimal tag value available to allow the
 * registration of nbtags data in a continuous range.
 *
 * Note that this function must be called exactly the same way on all nodes to
 * make sure the tags are identical from one node to another.
 *
 *******************************************************************************
 *
 * @param[in] nbtags
 *          The number of tags required to register the sparse matrix or right
 *          hand side.
 *
 *******************************************************************************
 *
 * @return V, the minimal tag value to use. The range [V:V+nbtags-1] is booked.
 *
 ********************************************************************************/
int64_t
pastix_starpu_tag_book( __attribute__((unused)) int64_t nbtags ) {
    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Release the set of tags starting by min.
 *
 * This function releases the range of tags that starts by the min value.
 *
 *******************************************************************************
 *
 * @param[in] min
 *          The initial value in the range
 *
 ******************************************************************************/
void
pastix_starpu_tag_release( __attribute__((unused)) int64_t min ) {
    return;
}

#endif

/**
 * @}
 */
