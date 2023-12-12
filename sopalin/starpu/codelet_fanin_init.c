/**
 *
 * @file codelet_fanin_init.c
 *
 * StarPU codelets to initialize a fanin before use. This codelet is not meant
 * to be submitted by the user, but automatically by the runtime at need.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.1
 * @author Alycia Lisito
 * @date 2023-12-01
 *
 * @addtogroup pastix_starpu
 * @{
 *
 **/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define _GNU_SOURCE
#include "common.h"
#include "pastix_starpu.h"
#include "codelets.h"

#if !defined(PASTIX_STARPU_SIMULATION)
/**
 *******************************************************************************
 *
 * @brief StarPU CPU implementation
 *
 *******************************************************************************
 *
 * @param[in] descr
 *          TODO
 *
 * @param[in] cl_arg
 *          TODO
 *
 *******************************************************************************/
static void
fct_fanin_init_cpu( void *descr[], void *cl_arg )
{
    void   *A    = pastix_starpu_cblk_get_ptr( descr[0] );
    size_t  size = ((pastix_starpu_interface_t *)descr[0])->allocsize;

    assert( size > 0 );
    memset( A, 0, size );

    (void)cl_arg;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

/**
 * @brief Main structure for all tasks of fanin_init type
 */
struct starpu_codelet cl_fanin_init_cpu = {
    .where     = STARPU_CPU,
#if !defined(PASTIX_STARPU_SIMULATION)
    .cpu_funcs[0] = fct_fanin_init_cpu,
#else
    .cpu_funcs[0] = (starpu_cpu_func_t)1,
#endif
    .nbuffers  = 1,
    .modes[0]  = STARPU_W,
    .name      = "fanin_init",
};

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @}
 */
