
/************************************************************/
/**                                                        **/
/**   NAME       : amalgamate.h                            **/
/**                                                        **/
/**   AUTHOR     : Mathieu Faverge                         **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 15/08/2006      **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#ifndef AMALGAMATE_H
#define AMALGAMATE_H

void amalgamate(double rat, csptr P, pastix_int_t snodenbr, pastix_int_t *snodetab, 
		pastix_int_t *treetab, pastix_int_t *cblknbr, pastix_int_t **rangtab, pastix_int_t *nodetab, MPI_Comm pastix_comm);

#endif
