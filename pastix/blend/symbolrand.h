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
PASTIX_INT hazard(PASTIX_INT a, PASTIX_INT b);
void symbolRand(SymbolMatrix *symbmtx, PASTIX_INT h1, PASTIX_INT h2);

#undef static
