/**
 *
 * @file pastix_task_symbfact.h
 *
 *  PaStiX symbolic factorizations routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains wrappers to the symbolic factorization step.
 * Affetcted by the compilation time options:
 *    - PASTIX_SYMBOL_FORCELOAD: Force to load the symbol matrix from file
 *    - PASTIX_SYMBOL_DUMP_SYMBMTX: Dump the symbol matrix in a postscript file.
 *    - COMPACT_SMX: Optimization for solve step (TODO: check if not obsolete)
 *    - FORGET_PARTITION: Force to forget the precomputed partition
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/

#ifndef PASTIX_TASK_SYMBFACT_H
#define PASTIX_TASK_SYMBFACT_H

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbfact
 * @ingroup pastix
 *
 * pastix_task_symbfact - Computes the the symbolic factorization of the matrix
 * and if required the amalgamated supernode partition.
 *
 * The function is a *centralized* algorithm to generate the symbol matrix
 * structure associated to the problem. It takes as input the ordemesh structure
 * (permutaion array, inverse permutation array, and optionnal supernodes
 * array) and returns the modified ordemesh structure if changed, and the
 * symbolic structure.
 *  - If (PT-)Scotch has been used, it generates the structure with
 * symbolFaxGraph() thanks to the supernode partition given by Scotch.
 *  - If ILU(k) factorization will be performed or if the ordering tools didn't
 * provide the supernode partition, symbolKass() is used to generate both
 * supernode partition and associated symbol matrix structure.
 *
 * Both algorithms are working with a centralized version of the graph and are
 * on every nodes. If a distributed graph has been used, it is gather on each
 * node to compute the symbol matrix.
 * If symbolKass() is used, the perm and invp vector will be modified and
 * returned to the user. BE CAREFULL if you give your own ordering and wants to
 * keep both version because the one given will be overwritten.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_INCOMPLETE, IPARM_LEVEL_OF_FILL,
 *   IPARM_AMALGAMATION_LVLCBLK, IPARM_AMALGAMATION_LVLBLAS,
 *   IPARM_IO_STRATEGY
 *
 *******************************************************************************
 *
 * @param[in,out] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the field symbmtx is initialized with the symbol matrix,
 *          and the field ordemesh is updated if the supernode partition is
 *          computed.
 *          - IPARM_INCOMPLETE switches the factorization mode from direct to ILU(k).
 *          - IPARM_LEVEL_OF_FILL defines the level of incomplete factorization
 *            if IPARM_INCOMPLETE == API_YES. If IPARM_LEVEL_OF_FILL < 0, the
 *            full pattern is generated as for direct factorization.
 *          - IPARM_AMALGAMATION_LVLCBLK is the ratio of amalgamation allowed
 *            based on reducing the number of supernodes only.
 *          - IPARM_AMALGAMATION_LVLBLAS is the ratio of amalgamation allowed
 *            based on reducing the computational cost (solve for ILU(k), or
 *            factorization for direct factorization).
 *          - IPARM_IO_STRATEGY will enable to load/store the result to files.
 *          If set to API_IO_SAVE, the symbmtx and the generated ordemesh is dump to file.
 *          If set to APÏ_IO_LOAD, the symbmtx (only) is loaded from the files.
 *
 * @param[in,out] perm
 *          Array of size n.
 *          On entry, unused.
 *          On exit, if perm != NULL, contains the permutation array generated.
 *
 * @param[in,out] invp
 *          Array of size n.
 *          On entry, unused.
 *          On exit, if invp != NULL, contains the inverse permutation array generated.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS on successful exit
 *          \retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *          \retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *          \retval PASTIX_ERR_INTEGER_TYPE if Scotch integer type is not the
 *                  same size as PaStiX ones.
 *          \retval PASTIX_ERR_INTERNAL if an error occurs internally to Scotch.
 *
 *******************************************************************************/
int
pastix_task_symbfact(d_pastix_data_t *pastix_data,
                     pastix_int_t  *perm,
                     pastix_int_t  *invp );
#endif /* PASTIX_TASK_SYMBFACT_H */
