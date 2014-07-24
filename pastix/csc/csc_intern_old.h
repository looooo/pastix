/*== Creation/Destruction de CSC ==*/

/******************************************************************************/
/* void CscHbRead(z_CscMatrix *thecsc, char *Type, char *RhsType, pastix_complex64_t **rhs,  */
/*                pastix_complex64_t **rhs2, const Order *ord, const char *rsaname)        */
/*                                                                            */
/* Lecture de la csc a partir d'un fichier au format HB et permutation a      */
/* partir du vecteur permutation fournit par Scotch                           */
/*                                                                            */
/* thecsc : La csc                                                            */
/* Type : type du HB (RUA, RSA ....)                                          */
/* RhsType : type pour les seconds membres                                    */
/* rhs : Vecteur second membre                                                */
/* rhs2 : Vecteur solution                                                    */
/* ord : la permutation                                                       */
/* rsaname : nom du fichier HB                                                */
/*                                                                            */
/* Type doit etre alloue avant l'appel, char Type[4]                          */
/* RhsType doit etre alloue avant l'appel, char RhsType[4]                    */
/* rhs et rhs2 sont alloue si besoin est.                                     */
/* RhsType[0] == '\0' si pas de second membre dans le fichier                 */
/******************************************************************************/

/* !!!!!!!!!!!!!FONCTION INUTILISEE !!!!!!!! */
void CscOrder(z_CscMatrix *thecsc,
	      char *Type, char *RhsType, pastix_complex64_t **transcsc,
	      pastix_complex64_t **rhs, pastix_complex64_t **rhs2,
	      const Order *ord,
	      pastix_int_t Nrow, pastix_int_t Ncol, pastix_int_t Nnzero, 
	      pastix_int_t *colptr, pastix_int_t *rowind, pastix_complex64_t *val, pastix_int_t forcetrans);

/*== Distribution/Remplissage ==*/
/******************************************************************************/
/* void CscDistrib(SymbolMatrix *symbmtx, z_CscMatrix *cscmtx,                  */
/*                 z_CscMatrix *cscdist)                                        */
/*                                                                            */
/* Distribution de la csc                                                     */
/*                                                                            */
/* symbmtx : symbol matrix locale                                             */
/* cscmtx : csc globale                                                       */
/* cscdist : csc locale                                                       */
/******************************************************************************/
void CscDistrib(const SymbolMatrix *symbmtx, const z_CscMatrix *cscmtx,
		z_CscMatrix *cscdist, const pastix_complex64_t *transcsc, pastix_complex64_t **trandcsc);

/*== Divers ==*/
/******************************************************************************/
/* void CscTrans(z_CscMatrix *cscmtx, z_CscMatrix *csctrp)                        */
/*                                                                            */
/* Transpose une csc                                                          */
/*                                                                            */
/* cscmtx : csc                                                               */
/* csctrp : csc transposee                                                    */
/******************************************************************************/
void CscTrans(const z_CscMatrix *cscmtx, z_CscMatrix *csctrp);

void z_CscScaling2(char *Type, pastix_int_t Ncol, pastix_int_t *col, pastix_int_t *row, pastix_complex64_t *val, pastix_complex64_t *rhs, pastix_complex64_t *rhs2);
void z_CscScaling(z_CscMatrix *cscmtx, pastix_complex64_t *transcsc, pastix_complex64_t *rhs1, pastix_complex64_t *rhs2);

/******************************************************************************/
/* void CscVerifUpdown(const d_UpDownVector *updovct,                           */
/*                     const SymbolMatrix *symbmtx; const pastix_complex64_t *rhs2)        */
/*                                                                            */
/* Verification entre le second membre fournit dans le fichier HB et le second*/
/* membre calcule.                                                            */
/*                                                                            */
/* updovct : vecteur second membre calcule                                    */
/* symbmtx : Symbol matrix                                                    */
/* rhs2 : vecteur second membre solution fournit dans le fichier HB           */
/******************************************************************************/
void CscVerifUpdown(const d_UpDownVector *updovct, const SymbolMatrix *symbmtx,
		    const pastix_complex64_t *rhs2);

/******************************************************************************/
/* void CscUpDown(d_UpDownVector *updovct, const SymbolMatrix *symbmtx,         */
/*                const pastix_complex64_t *rhs)                                           */
/*                                                                            */
/* Remplissage du vector second membre a partir de celui fournit dans le      */
/* fichier HB.                                                                */
/*                                                                            */
/******************************************************************************/
void CscUpdown(d_UpDownVector *updovct, /*const*/ SymbolMatrix *symbmtx,
	       const pastix_complex64_t *rhs);
void CscUpdown2(d_UpDownVector *updovct, /*const*/ SymbolMatrix *symbmtx,
		const pastix_complex64_t *rhs);



/******************************************************************************/
/* void CscDiagDom(z_CscMatrix *cscmtx)                                         */
/*                                                                            */
/* Transforme la csc en csc a diagonale dominante                             */
/*                                                                            */
/* cscmtx : Csc                                                               */
/******************************************************************************/
void CscDiagDom(z_CscMatrix *cscmtx);

pastix_int_t CscStrucSym(z_CscMatrix *cscmtx);
