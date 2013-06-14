
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

void amalgamate(double rat, csptr P, PASTIX_INT snodenbr, PASTIX_INT *snodetab, 
		PASTIX_INT *treetab, PASTIX_INT *cblknbr, PASTIX_INT **rangtab, PASTIX_INT *nodetab, MPI_Comm pastix_comm);

#endif
