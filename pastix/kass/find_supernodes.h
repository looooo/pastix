/************************************************************/
/**                                                        **/
/**   NAME       : find_supernodes.h                       **/
/**                                                        **/
/**   AUTHOR     : Mathieu Faverge                         **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10/02/2006      **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#ifndef FIND_SUPERNODES_H
#define FIND_SUPERNODES_H

void  find_supernodes(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_int_t *perm, pastix_int_t *iperm, 
		      pastix_int_t *snodenbr, pastix_int_t *snodetab, pastix_int_t *treetab);

#endif


