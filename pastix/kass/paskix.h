/************************************************************/
/**                                                        **/
/**   NAME       :  paskix.h                              **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 28 Nov 2006     **/
/**                                                        **/
/**  Some functions on CSC matrices                        **/
/**                                                        **/
/************************************************************/

/*
**  The function prototypes.
*/
void CSC_Fnum2Cnum(PASTIX_INT *ja, PASTIX_INT *ia, PASTIX_INT n);
void CSC_Cnum2Fnum(PASTIX_INT *ja, PASTIX_INT *ia, PASTIX_INT n);

void CSC_rowPerm(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a, PASTIX_INT *rperm);
void CSC_colPerm(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a, PASTIX_INT *cperm);
void CSC_colScale(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a, PASTIX_FLOAT *dcol);
void CSC_rowScale(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a, PASTIX_FLOAT *drow);
void CSC_sort(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a);
