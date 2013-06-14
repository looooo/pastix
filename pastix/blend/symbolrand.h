/************************************************************/
/**                                                        **/
/**   NAME       : symbolrand.h                            **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : ILU incomplet a la louche               **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#define static
pastix_int_t hazard(pastix_int_t a, pastix_int_t b);
void symbolRand(SymbolMatrix *symbmtx, pastix_int_t h1, pastix_int_t h2);

#undef static
