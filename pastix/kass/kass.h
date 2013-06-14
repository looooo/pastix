/************************************************************/
/**                                                        **/
/**   NAME       : kass.h                                  **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Compute a block structure of the factor **/
/**                obtained by a ILU(k) factorization      **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 30/01/2006      **/
/**                                 to                     **/
/**                                                        **/
/************************************************************/

void kass(int            levelk, 
	  int            rat, 
	  SymbolMatrix * symbptr, 
	  PASTIX_INT            baseval,
	  PASTIX_INT            vertnbr, 
	  PASTIX_INT            edgenbr, 
	  PASTIX_INT          * verttab,
	  PASTIX_INT          * edgetab, 
	  Order        * orderptr, 
	  MPI_Comm       pastix_comm);

/* void kass(int alpha, int rat, SymbolMatrix * symbptr, Graph * graphptr, Order * orderptr, MPI_Comm pastix_comm); */

void ifax(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_INT levelk, PASTIX_INT  cblknbr, PASTIX_INT *rangtab, PASTIX_INT *perm, PASTIX_INT *iperm, SymbolMatrix *symbmtx);

