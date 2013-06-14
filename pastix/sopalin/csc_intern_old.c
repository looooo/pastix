/*== Creation/Destruction de CSC ==*/

/******************************************************************************/
/* void CscOrder(CscMatrix *thecsc, char *Type, char *RhsType, pastix_float_t **rhs,  */
/*                pastix_float_t **rhs2, const Order *ord, ...)                        */
/*                                                                            */
/* Construction de la csc a partir de col/row/val et permutation a            */
/* partir du vecteur permutation fournit par Scotch                           */
/*                                                                            */
/* thecsc : La csc                                                            */
/* Type : type du HB (RUA, RSA ....)                                          */
/* RhsType : type pour les seconds membres                                    */
/* rhs : Vecteur second membre                                                */
/* rhs2 : Vecteur solution                                                    */
/* ord : la permutation                                                       */
/*                                                                            */
/* Type doit etre alloue avant l'appel, char Type[4]                          */
/* RhsType doit etre alloue avant l'appel, char RhsType[4]                    */
/* rhs et rhs2 sont alloue si besoin est.                                     */
/* RhsType[0] == '\0' si pas de second membre dans le fichier                 */
/******************************************************************************/

/* !!! FONCTION INUTILISEE !!! */
void CscOrder(CscMatrix *thecsc,
	      char *Type, char *RhsType, pastix_float_t **transcsc,
	      pastix_float_t **rhs, pastix_float_t **rhs2,
	      const Order *ord, 
	      pastix_int_t Nrow, pastix_int_t Ncol, pastix_int_t Nnzero, 
	      pastix_int_t *colptr, pastix_int_t *rowind, pastix_float_t *val, pastix_int_t forcetrans)
{
  char *crhs=NULL; /* 2nd member vector */
  char *crhs2=NULL; /* solution vector */ 

  pastix_int_t index,itercol,iter,newcol,colidx, rowp1; 
  pastix_int_t valnbr=0;
  pastix_int_t *trscltb=NULL; /* coltab for transpose csc */
  pastix_int_t *trowtab=NULL; /* rowtab for transpose csc */

#ifdef CSC_LOG
  fprintf(stdout, "-> CscOrder \n");
#endif

  /* Csc alloc */
  CSC_FNBR(thecsc) = 1;
  CSC_FTAB(thecsc) = (CscFormat *) memAlloc(1*sizeof(CscFormat));
  if (CSC_FTAB(thecsc) == NULL)
    errorPrint( "CscHbRead : Not enough memory for CSC_FTAB\n");
  CSC_COLNBR(thecsc,0) = Ncol;
  CSC_COLTAB(thecsc,0) = (pastix_int_t *) memAlloc((Ncol+1)*sizeof(pastix_int_t));
  if (CSC_COLTAB(thecsc,0) == NULL)
    errorPrint( "CscHbRead : Not enough memory for CSC_COLTAB\n");

  for (index=0; index < (Ncol+1); index++)
    CSC_COL(thecsc,0,index) = 0;

  if (Type[1] == 'S')
    {
      /* Symetric */
      CSC_ROWTAB(thecsc) = (pastix_int_t*) memAlloc(2*Nnzero*sizeof(pastix_int_t));
      if (CSC_ROWTAB(thecsc) == NULL)
	errorPrint( "CscHbRead : Not enough memory for CSC_ROWTAB\n");
      CSC_VALTAB(thecsc) = (pastix_float_t*) memAlloc(2*Nnzero*sizeof(pastix_float_t));
      if (CSC_VALTAB(thecsc) == NULL)
	errorPrint( "CscHbRead : Not enough memory for CSC_VALTAB\n");
      if (forcetrans)
	{
	  printf("Force transpose on symetric...\n");

	  (*transcsc) = (pastix_float_t *) memAlloc(2*Nnzero*sizeof(pastix_float_t));
	  if ((*transcsc) == NULL)
	    errorPrint( "CscHbRead : Not enough memory for (*transcsc)\n");
	  trowtab = (pastix_int_t *) memAlloc(2*Nnzero*sizeof(pastix_int_t));
	  if (trowtab == NULL)
	    errorPrint( "CscHbRead : Not enough memory for trowtab\n");
	}
    }
  else
    {
      /* Unsymmetric */
      CSC_ROWTAB(thecsc) = (pastix_int_t*) memAlloc(Nnzero*sizeof(pastix_int_t));
      if (CSC_ROWTAB(thecsc) == NULL)
	errorPrint( "CscHbRead : Not enough memory for CSC_ROWTAB\n");
      CSC_VALTAB(thecsc) = (pastix_float_t*) memAlloc(Nnzero*sizeof(pastix_float_t));
      if (CSC_VALTAB(thecsc) == NULL)
	errorPrint( "CscHbRead : Not enough memory for CSC_VALTAB\n");
      
      (*transcsc) = (pastix_float_t *) memAlloc(Nnzero*sizeof(pastix_float_t));
      if ((*transcsc) == NULL)
	errorPrint( "CscHbRead : Not enough memory for (*transcsc)\n");
      trowtab = (pastix_int_t *) memAlloc(Nnzero*sizeof(pastix_int_t));
      if (trowtab == NULL)
	errorPrint( "CscHbRead : Not enough memory for trowtab\n");
    }

  /* Computing good coltabs */
  for (itercol = 0; itercol < Ncol; itercol++)
    {
      newcol = ord->permtab[itercol];

      CSC_COL(thecsc,0,newcol) += colptr[itercol+1] - colptr[itercol];

      if (Type[1] == 'S')
	{
	  /* Symmetric */
	  for (iter=colptr[itercol]; iter<colptr[itercol+1]; iter++)
	    {
	      if ((rowind[iter-1]-1) != itercol)
		{
		  newcol = ord->permtab[rowind[iter-1]-1];
		  (CSC_COL(thecsc,0,newcol))++;
		}
	    }
	}
    }

  newcol = 0;
  for (index=0; index<(Ncol+1); index++)
    {
      colidx = CSC_COL(thecsc,0,index);
      CSC_COL(thecsc,0,index) = newcol;
      newcol += colidx;
    }

  if ((*transcsc) != NULL)
    {
      trscltb = (pastix_int_t *) memAlloc((Ncol+1)*sizeof(pastix_int_t));
      if (trscltb == NULL)
	errorPrint( "CscHbRead : not enough memory for trscltb\n");
      for (index=0; index<(Ncol+1); index++)
	{
	  trscltb[index] = CSC_COL(thecsc,0,index);
	}
    }

  /* Put the element */

  printf("Debut construction CSC\n");
  for (itercol = 0; itercol <Ncol; itercol++)
    {
      for (iter=colptr[itercol]; iter<colptr[itercol+1]; iter++)
	{
	  pastix_int_t therow;

	  rowp1 = rowind[iter-1] -1;
	  newcol = ord->permtab[itercol];
	  colidx = CSC_COL(thecsc,0,newcol);
	  therow = ord->permtab[rowp1];
	  CSC_ROW(thecsc,colidx) = therow;
	  CSC_VAL(thecsc,colidx) = val[iter-1];	  
	  valnbr++;

	  CSC_COL(thecsc,0,newcol)++;

	  if ((*transcsc) != NULL)
	    {
	      (*transcsc)[trscltb[therow]] = val[iter-1];	      
	      trowtab[trscltb[therow]] = newcol;
	      trscltb[therow]++;
	    }


	  /* Symmetric */
	  if (Type[1] == 'S')
	    {
	      if (rowp1 != itercol)
		{
		  newcol = ord->permtab[rowp1];
		  colidx = CSC_COL(thecsc,0,newcol);
		  therow = ord->permtab[itercol];
		  CSC_ROW(thecsc,colidx) = therow;
		  CSC_VAL(thecsc,colidx) = val[iter-1];
		  valnbr++;

		  (CSC_COL(thecsc,0,newcol))++;
		  
		  if ((*transcsc) != NULL)
		    {
		      (*transcsc)[trscltb[therow]] = val[iter-1];
		      trowtab[trscltb[therow]] = newcol;
		      trscltb[therow]++;
		    }
		}
	    }
	  
	}
    }
  memFree_null(colptr);
  memFree_null(rowind);
  memFree_null(val);
  if (trscltb != NULL)
    memFree_null(trscltb);
  printf("Fin construction CSC\n");

  /* 2nd member */
  if (RhsType[0] != '\0')
    {
      printf("Not yet implemented\n");
      EXIT(MOD_SI,NOTIMPLEMENTED_ERR);

      (*rhs) = (pastix_float_t*) memAlloc(Ncol*sizeof(pastix_float_t));
      if ((*rhs) == NULL)
	errorPrint( "CscHbRead : Not enough memory for rhs\n");
      
      for (index=0; index<Ncol; index++)
	{
	  (*rhs)[ord->permtab[index]] = crhs[index];
	}

      if (RhsType[2] == 'X')
	{
	  /* Vector Solution */
	  (*rhs2) = (pastix_float_t*) memAlloc(Ncol*sizeof(pastix_float_t));
	  if ((*rhs2) == NULL)
	    errorPrint( "CscHbRead : not enough memory for rhs2\n");

	  for (index=0; index<Ncol; index++)
	    {
	      (*rhs2)[ord->permtab[index]] = crhs2[index];
	    }
	}
      
      memFree_null(crhs);
      memFree_null(crhs2);
    }
  else
    {
      (*rhs) = NULL;
      (*rhs2) = NULL;
    }
  printf("valnbr = %ld\n", (long)valnbr);

  /* good coltab */
  colidx = 0;
  for (index=0; index<Ncol; index++)
    {
      newcol = CSC_COL(thecsc,0,index);
      CSC_COL(thecsc,0,index) = colidx;
      colidx = newcol;
    }

  /* Sort on the row */
  for (index = 0; index < Ncol; index++)
    {
      /* bubble sort between coltab[index2] and coltab[index2+1] */
      pastix_int_t *t = &(CSC_FROW(thecsc,0,index));
      pastix_float_t *v = &(CSC_FVAL(thecsc,0,index));
      pastix_int_t n = CSC_COL(thecsc,0,index+1) - CSC_COL(thecsc,0,index);
	
      pastix_int_t i,j;
      for (i=0; i<n; i++)
	for (j=0; j<n-i-1; j++)
	  if (t[j] > t[j+1])
	    {
	      pastix_int_t tempt = t[j+1];
	      pastix_float_t tempv = v[j+1];
	      
	      t[j+1] = t[j];
	      v[j+1] = v[j];
	      t[j] = tempt;
	      v[j] = tempv;
	    }
    }

  if ((*transcsc) != NULL)
    {
      for (index=0; index<Ncol; index++)
	{
	  pastix_int_t *t = &(trowtab[CSC_COL(thecsc,0,index)]);
	  pastix_float_t *v = &((*transcsc)[CSC_COL(thecsc,0,index)]);
	  
	  pastix_int_t n = CSC_COL(thecsc,0,index+1) - CSC_COL(thecsc,0,index);
	  pastix_int_t i,j;
	  
	  for (i=0; i<n; i++)
	    for (j=0; j<n-i-1; j++)
	      if (t[j] > t[j+1])
		{
		  pastix_int_t tempt = t[j+1];
		  pastix_float_t tempv = v[j+1];
		  
		  t[j+1] = t[j];
		  v[j+1] = v[j];
		  t[j] = tempt;
		  v[j] = tempv;
		}
	}
      memFree_null(trowtab);
    }
#ifdef CSC_LOG
  fprintf(stdout, "<- CscOrder \n");
#endif
}

/*== Distribution/Remplissage ==*/
/******************************************************************************/
/* void CscDistrib(SymbolMatrix *symbmtx, CscMatrix *cscmtx,                  */
/*                 CscMatrix *cscdist, const pastix_float_t *transcsc,                 */
/*		   pastix_float_t **trandcsc)                                          */
/*                                                                            */
/* Distribution de la csc                                                     */
/*                                                                            */
/* symbmtx : symbol matrix locale                                             */
/* cscmtx : csc globale                                                       */
/* cscdist : csc locale                                                       */
/* transcsc : csr globale                                                     */
/* trandcsc : csr locale                                                      */
/******************************************************************************/
/* !!! FONCTION INUTILISEE !!! */
void CscDistrib(const SymbolMatrix *symbmtx, const CscMatrix *cscmtx,
		CscMatrix * cscdist, const pastix_float_t *transcsc, pastix_float_t **trandcsc)
{
  pastix_int_t iterckcd=0; /* iter for cblk in csc(d) */
  pastix_int_t itercol; /* iter for coltab */
  pastix_int_t strdcol=0; /* stride for col */
  pastix_int_t itervald=0; /* iter for valtab rowtab */
  pastix_int_t itervalg=0;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscDistrib \n");
#endif
  
  CSC_FNBR(cscdist) = symbmtx->cblknbr;
  CSC_FTAB(cscdist) = (CscFormat *) memAlloc(CSC_FNBR(cscdist)*sizeof(CscFormat));
  if (CSC_FTAB(cscdist) == NULL)
    errorPrint( "CscDistrib : Not enough memory for CSC_FTAB\n");
  
  
  for (iterckcd=0;
       iterckcd < symbmtx->cblknbr;
       iterckcd++)
    {
      CSC_COLNBR(cscdist,iterckcd) =
	symbmtx->cblktab[iterckcd].lcolnum-symbmtx->cblktab[iterckcd].fcolnum+1;
      CSC_COLTAB(cscdist,iterckcd) =
	(pastix_int_t *) memAlloc((CSC_COLNBR(cscdist,iterckcd)+1)*sizeof(pastix_int_t));
      if (CSC_COLTAB(cscdist,iterckcd) == NULL)
	errorPrint( "CscDistrib : Not enough memory for CSC_COLTAB\n");

      for (itercol = 0;
	   itercol < (CSC_COLNBR(cscdist,iterckcd)+1);
	   itercol++)
	{
	  CSC_COL(cscdist,iterckcd,itercol) =
	    CSC_COL(cscmtx,0,itercol+symbmtx->cblktab[iterckcd].fcolnum) -
	    CSC_COL(cscmtx,0,symbmtx->cblktab[iterckcd].fcolnum) +
	    strdcol;
	}
      strdcol = CSC_COL(cscdist,iterckcd,CSC_COLNBR(cscdist,iterckcd));
    }

  /* Remplissage des valtabs et rowtabs */
  CSC_ROWTAB(cscdist) = (pastix_int_t *) memAlloc(strdcol*sizeof(pastix_int_t));
  if (CSC_ROWTAB(cscdist) == NULL)
    errorPrint( "CscDistrib : Not enough memory for CSC_ROWTAB\n");
  CSC_VALTAB(cscdist) = (pastix_float_t *) memAlloc(strdcol*sizeof(pastix_float_t));
  if (CSC_VALTAB(cscdist) == NULL)
    errorPrint( "CscDistrib : Not enough memory for CSC_VALTAB\n");

  if (transcsc != NULL)
    {
      /* Rua */
      (*trandcsc) = (pastix_float_t *) memAlloc(strdcol*sizeof(pastix_float_t));
      if ((*trandcsc) == NULL)
	errorPrint( "CscDistrib : Not enough memory for trandcsc\n");
    }
  
  for (iterckcd = 0;
       iterckcd < CSC_FNBR(cscdist);
       iterckcd++)
    {
      for (itercol = 0;
	   itercol < CSC_COLNBR(cscdist,iterckcd);
	   itercol++)
	{
	  itervalg =
	    CSC_COL(cscmtx,0,itercol+symbmtx->cblktab[iterckcd].fcolnum);
	  
	  for (itervald = CSC_COL(cscdist,iterckcd,itercol);
	       itervald < CSC_COL(cscdist,iterckcd,itercol+1);
	       itervald++)
	    {
	      CSC_VAL(cscdist,itervald) = CSC_VAL(cscmtx,itervalg);

	      if (transcsc != NULL)
		(*trandcsc)[itervald] = transcsc[itervalg];
	      
	      CSC_ROW(cscdist,itervald) = CSC_ROW(cscmtx,itervalg);
	      itervalg++;
	    }
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscDistrib \n");
#endif
}


/******************************************************************************/
/* void Csc2solv(CscMatrix *cscmtx, SolverMatrix *solvmtx, pastix_float_t *trandcsc)   */
/*                                                                            */
/* Remplit la solvermatrix locale a partir de la csc locale, pour la partie   */
/* on utiliser la transpose de la csc globale qui est lu a partir d'un fichier*/
/*                                                                            */
/* cscmtx : csc locale                                                        */
/* solvmtx : solver locale                                                    */
/* stream : fichier de la csc transposee ouvert en lecture                    */
/******************************************************************************/
/* !!! FONCTION INUTILISEE !!! */
void Csc2symb(const CscMatrix *cscmtx, SymbolMatrix *symbmtx)
{
  pastix_int_t itercblk;
  pastix_int_t itercoltab;
  pastix_int_t iterbloc;
  pastix_int_t iterval;

#ifdef CSC_LOG
  fprintf(stdout, "-> Csc2symb \n");
#endif
  
  for (iterbloc=0; iterbloc < symbmtx->bloknbr; iterbloc++)
    symbmtx->bloktab[iterbloc].levfval=0;

  for (itercblk=0; itercblk < CSC_FNBR(cscmtx); itercblk++)
    {
      for (itercoltab=0;
	   itercoltab < CSC_COLNBR(cscmtx,itercblk);
	   itercoltab++)
	{
	  for (iterval = CSC_COL(cscmtx,itercblk,itercoltab);
	       iterval < CSC_COL(cscmtx,itercblk,itercoltab+1);
	       iterval++)
	    {
	      if (CSC_ROW(cscmtx,iterval) >=
		  symbmtx->cblktab[itercblk].fcolnum)
		{
		  iterbloc = symbmtx->cblktab[itercblk].bloknum;

		  while ((( symbmtx->bloktab[iterbloc].lrownum <
			    CSC_ROW(cscmtx,iterval)) ||
			  ( symbmtx->bloktab[iterbloc].frownum >
			    CSC_ROW(cscmtx,iterval))) &&
			 ( iterbloc < symbmtx->cblktab[itercblk+1].bloknum))
		    {
		      iterbloc++;
		    }

		  if ( iterbloc <
		       symbmtx->cblktab[itercblk+1].bloknum)
		    {
		      symbmtx->bloktab[iterbloc].levfval=1;
		    }
		  else printf("ILU: csc2symb drop coeff from CSC c=%ld(%ld) l=%ld(%ld) cblk=%ld fcol=%ld lcol=%ld\n",
			      (long)symbmtx->cblktab[itercblk].fcolnum+
			      (long)itercoltab,(long)itercoltab,
			      (long)CSC_ROW(cscmtx,iterval),(long)iterval,
			      (long)itercblk,
			      (long)symbmtx->cblktab[itercblk].fcolnum,
			      (long)symbmtx->cblktab[itercblk].lcolnum);
		}
	    }
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- Csc2symb \n");
#endif
}


/*== Divers ==*/
/******************************************************************************/
/* void CscTrans(CscMatrix *cscmtx, CscMatrix *csctrp)                        */
/*                                                                            */
/* Transpose une csc                                                          */
/*                                                                            */
/* cscmtx : csc                                                               */
/* csctrp : csc transposee                                                    */
/******************************************************************************/
/* !!! FONCTION INUTILISEE !!! */
void CscTrans(const CscMatrix *cscmtx, CscMatrix *csctrp)
{
  pastix_int_t itercscf;
  pastix_int_t itercol;
  pastix_int_t valnbr;
  pastix_int_t colcur;
  pastix_int_t iterval;
  pastix_int_t rowcur;
  pastix_int_t iterval2;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscTrans \n");
#endif

  /* Copie du coltab */
  CSC_FNBR(csctrp) = CSC_FNBR(cscmtx);
  CSC_FTAB(csctrp) = (CscFormat *) memAlloc(CSC_FNBR(cscmtx)*sizeof(CscFormat));
  if (CSC_FTAB(csctrp) == NULL)
    errorPrint( "CscTrans : Not enough memory for CSC_FTAB\n");

  for (itercscf = 0;
       itercscf < CSC_FNBR(csctrp);
       itercscf++)
    {
      CSC_COLNBR(csctrp,itercscf) = CSC_COLNBR(cscmtx,itercscf);
      CSC_COLTAB(csctrp,itercscf) =
	(pastix_int_t *) memAlloc((CSC_COLNBR(csctrp,itercscf)+1)*sizeof(pastix_int_t));
      if (CSC_COLTAB(csctrp,itercscf) == NULL)
	errorPrint( "CscTrans : Not enough memory for CSC_COLTAB\n");
      
      for (itercol = 0;
	   itercol < (CSC_COLNBR(csctrp,itercscf)+1);
	   itercol++)
	{
	  CSC_COL(csctrp,itercscf,itercol) = CSC_COL(cscmtx,itercscf,itercol);
	}
    }

  /* Allocation du valtab et rowtab */
  valnbr = CSC_VALNBR(cscmtx);
  printf("valnbr = %ld\n", (long)valnbr);
  CSC_ROWTAB(csctrp) = (pastix_int_t *) memAlloc(valnbr*sizeof(pastix_int_t));
  if (CSC_ROWTAB(csctrp) == NULL)
    errorPrint( "CscTrans : Not enough memory for CSC_ROWTAB\n");
  CSC_VALTAB(csctrp) = (pastix_float_t *) memAlloc(valnbr*sizeof(pastix_float_t));
  if (CSC_VALTAB(csctrp) == NULL)
    errorPrint( "CscTrans : Not enough memory for CSC_VALTAB\n");

  /* Renseignement des coeff */
  colcur = 0;
  for (itercscf = 0;
       itercscf < CSC_FNBR(cscmtx);
       itercscf++)
    {
      for (itercol = 0;
	   itercol < CSC_COLNBR(cscmtx,itercscf);
	   itercol++)
	{
	  for (iterval = CSC_COL(cscmtx,itercscf,itercol);
	       iterval < CSC_COL(cscmtx,itercscf,itercol+1);
	       iterval++)
	    {
	      pastix_int_t cont = 0;
	      pastix_int_t itercscf2 = 0;

	      rowcur = CSC_ROW(cscmtx,iterval);

	      /* Recherche dur coltab correspondant a rowcur dans trp */
	      do
		{
		  if (rowcur < CSC_COLNBR(csctrp,itercscf2))
		    {
		      iterval2 = CSC_COL(csctrp,itercscf2,rowcur);

		      CSC_ROW(csctrp,iterval2) = colcur;

		      CSC_VAL(csctrp,iterval2) = CSC_VAL(cscmtx,iterval);

		      (CSC_COL(csctrp,itercscf2,rowcur))++;

		      cont = 1;
		    }
		  else
		    {
		      rowcur -= CSC_COLNBR(csctrp,itercscf2);
		      itercscf2++;
		    }
		}
	      while (cont==0);
	    }
	  colcur++;
	}
    }

  rowcur = 0;
  /* Remet le coltab */
  for (itercscf = 0;
       itercscf < CSC_FNBR(csctrp);
       itercscf++)
    {
      for (itercol = 0;
	   itercol < (CSC_COLNBR(csctrp,itercscf));
	   itercol++)
	{
	  pastix_int_t rowidx = CSC_COL(csctrp,itercscf,itercol);
	  CSC_COL(csctrp,itercscf,itercol) = rowcur;
	  rowcur = rowidx;
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscTrans \n");
#endif
}

/* !!! FONCTION INUTILISEE !!! */
void CscScaling2(char *Type, pastix_int_t Ncol, pastix_int_t *col, pastix_int_t *row, pastix_float_t *val, pastix_float_t *rhs, pastix_float_t *rhs2)
{
  pastix_int_t itercol;
  pastix_int_t iterval;
  pastix_float_t *rowscal;
  pastix_float_t *colscal;
  const pastix_int_t flag = (Type[1] != 'S');

#ifdef CSC_LOG
  fprintf(stdout, "-> CscScaling2 \n");
#endif

  if (flag)
    rowscal = (pastix_float_t *) memAlloc(Ncol*sizeof(pastix_float_t));
  else
    rowscal = NULL;
  colscal = (pastix_float_t *) memAlloc(Ncol*sizeof(pastix_float_t));

  if (flag)
    if (rowscal == NULL)
      errorPrint( "CscScaling2 : Not enough memory for rowscal\n");
  if (colscal == NULL)
    errorPrint( "CscScalign2 : Nto enough memory for colscal\n");
  
  for (itercol=0; itercol<Ncol; itercol++)
    {
      if (flag)
	rowscal[itercol] = 0;
      colscal[itercol] = 0;
    }

  for (itercol=0; itercol<Ncol; itercol++)
    {
      for (iterval=col[itercol]; iterval<col[itercol+1]; iterval++)
	{
	  colscal[itercol] += ABS_FLOAT(val[iterval-1]);
	  if (flag)
	    rowscal[row[iterval-1]-1] += ABS_FLOAT(val[iterval-1]);
	}
    }

  for (itercol=0; itercol<Ncol; itercol++)
    {
      for (iterval=col[itercol]; iterval<col[itercol+1]; iterval++)
	{
	  val[iterval-1] /= colscal[itercol];
	  if (flag)
	    val[iterval-1] /= rowscal[row[iterval-1]-1];
	  else
	    val[iterval-1] /= colscal[itercol];
	}
      if (rhs != NULL)
	{
	  if (flag)
	    rhs[itercol] /= rowscal[itercol];
	  else
	    rhs[itercol] /= colscal[itercol];
	}
      if (rhs2 != NULL)
	rhs2[itercol] *= colscal[itercol];
    }

  if (flag)
    memFree_null(rowscal);
  memFree_null(colscal);

#ifdef CSC_LOG
  fprintf(stdout, "<- CscScaling2 \n");
#endif
}

/******************************************************************************/
/* void CscScaling(CscMatrix *cscmtx, pastix_float_t *transcsc,                        */
/*                 pastix_float_t *rhs, pastix_float_t *rhs2)                                   */
/*                                                                            */
/* Scaling                                                                    */
/*                                                                            */
/* cscmtx : csc                                                               */
/* transcsc : transpose csc (NULL si on est en RSA)                           */
/* rhs : second membre                                                        */
/* rhs2 : solution                                                            */
/******************************************************************************/
/* !!! FONCTION INUTILISEE !!! */
void CscScaling(CscMatrix *cscmtx, pastix_float_t *transcsc, pastix_float_t *rhs, pastix_float_t *rhs2)
{
  const pastix_int_t itercscf=0;
  pastix_int_t itercol;
  pastix_int_t iterval;
  pastix_float_t *rowscal;
  pastix_float_t *colscal;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscScaling \n");
#endif
  
  rowscal = (pastix_float_t *) memAlloc(CSC_COLNBR(cscmtx,0)*sizeof(pastix_float_t));
  colscal = (pastix_float_t *) memAlloc(CSC_COLNBR(cscmtx,0)*sizeof(pastix_float_t));
  
  for (itercol=0; itercol<CSC_COLNBR(cscmtx,itercscf); itercol++)
    {
      rowscal[itercol] = 0;
      colscal[itercol] = 0;
    }

  /* Compute Scaling  */
  for (itercol=0; itercol<CSC_COLNBR(cscmtx,itercscf); itercol++)
    {
      for (iterval=CSC_COL(cscmtx,itercscf,itercol);
	   iterval<CSC_COL(cscmtx,itercscf,itercol+1);
	   iterval++)
	{
	  colscal[itercol] += ABS_FLOAT(CSC_VAL(cscmtx,iterval));
	  rowscal[CSC_ROW(cscmtx,iterval)] += ABS_FLOAT(CSC_VAL(cscmtx,iterval));
	}
    }

  /* Apply Scaling */
  for (itercol=0; itercol<CSC_COLNBR(cscmtx,itercscf); itercol++)
    {
      for (iterval=CSC_COL(cscmtx,itercscf,itercol);
	   iterval<CSC_COL(cscmtx,itercscf,itercol+1);
	   iterval++)
	{
	  CSC_VAL(cscmtx,iterval) /= colscal[itercol];
	  CSC_VAL(cscmtx,iterval) /= rowscal[CSC_ROW(cscmtx,iterval)];

	  if (transcsc != NULL)
	    {
	      transcsc[iterval] /= rowscal[itercol];
	      transcsc[iterval] /= colscal[CSC_ROW(cscmtx,iterval)];
	    }
	}
      if (rhs != NULL) rhs[itercol] /= rowscal[itercol];
      if (rhs2 != NULL) rhs2[itercol] *= colscal[itercol];
    }

  memFree_null(rowscal);
  memFree_null(colscal);

#ifdef CSC_LOG
  fprintf(stdout, "<- CscScaling \n");
#endif
}

/* !!! FONCTION INUTILISEE !!! */
void CscVerifUpdown(const UpDownVector *updovct, const SymbolMatrix *symbmtx,
		    const pastix_float_t *rhs2)
{
  pastix_int_t itercblk;
  pastix_int_t itercol;
  pastix_int_t itersm2x;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscVerifUpdown \n");
#endif

  for (itercblk=0; itercblk<symbmtx->cblknbr; itercblk++)
    {
      itersm2x = updovct->cblktab[itercblk].sm2xind;

      for (itercol=symbmtx->cblktab[itercblk].fcolnum;
	   itercol<symbmtx->cblktab[itercblk].lcolnum+1;
	   itercol++)
	{
#ifdef CPLX
	  errorPrint( "%ld (%10e,%10e) (%10e,%10e)\n",
		  (long)itercol, creal(rhs2[itercol]), cimag(rhs2[itercol]), creal(updovct->sm2xtab[itersm2x]), cimag(updovct->sm2xtab[itersm2x]));
#else
	  errorPrint( "%ld %10e %10e\n",
		  (long)itercol, rhs2[itercol], updovct->sm2xtab[itersm2x]);
#endif
	  itersm2x++;
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscVerifUpdown \n");
#endif
}

/* !!! FONCTION INUTILISEE !!! */
void CscUpdown(UpDownVector *updovct, /*const*/ SymbolMatrix *symbmtx,
	       const pastix_float_t *rhs)
{
  pastix_int_t itercblk;
  pastix_int_t itercol;
  pastix_int_t itersm2x;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscUpdown \n");
#endif
  
  for (itercblk=0; itercblk<symbmtx->cblknbr; itercblk++)
    {
      itersm2x = updovct->cblktab[itercblk].sm2xind;
      for (itercol=symbmtx->cblktab[itercblk].fcolnum;
	   itercol<symbmtx->cblktab[itercblk].lcolnum+1;
	   itercol++)
	{
	  updovct->sm2xtab[itersm2x] = rhs[itercol];
	  itersm2x++;
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscUpdown \n");
#endif
}

/* !!! FONCTION INUTILISEE !!! */
void CscUpdown2(UpDownVector *updovct, /*const*/ SymbolMatrix *symbmtx,
		const pastix_float_t *rhs)
{
  pastix_int_t itercblk;
  pastix_int_t itercol;
  pastix_int_t itersm2x;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscUpdown2 \n");
#endif

  for (itercblk=0; itercblk<symbmtx->cblknbr; itercblk++)
    {
      itersm2x = updovct->cblktab[itercblk].sm2xind;
      for (itercol=symbmtx->cblktab[itercblk].fcolnum;
	   itercol<symbmtx->cblktab[itercblk].lcolnum+1;
	   itercol++)
	{
	  updovct->sm2xtab[itersm2x] = 1.0;
	  itersm2x++;
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscUpdown2 \n");
#endif
}
/* !!! INUTILISEE, non testee !!! */
void Csc2updown_new(Sopalin_Data_t * sopalin_data, int me, const CscMatrix *cscmtx, UpDownVector *updovct,
		    /*const*/ SymbolMatrix *symbmtx, int n, MPI_Comm comm)
{
  SolverMatrix * datacode;
  pastix_float_t *tempy;
  pastix_int_t    i;

#ifdef CSC_LOG
  fprintf(stdout, "-> Csc2updown_new \n");
#endif

  datacode = sopalin_data->datacode;

  MONOTHREAD_BEGIN;

  if (!(tempy = (pastix_float_t *)memAlloc(updovct->gnodenbr*sizeof(pastix_float_t)))) MALLOC_ERROR("tempy");
  sopalin_data->ptr_raff[0] = (void *)tempy;

  for (i=0;i<updovct->gnodenbr;i++)
    tempy[i] = (pastix_float_t)i;
  MONOTHREAD_END;

  SYNCHRO_THREAD;

  tempy = (pastix_float_t *)sopalin_data->ptr_raff[0];

  CscAx(sopalin_data, me, cscmtx, tempy, updovct->sm2xtab, symbmtx, updovct, comm);

  MONOTHREAD_BEGIN;
  memFree_null(tempy);
  MONOTHREAD_END;

#ifdef CSC_LOG
  fprintf(stdout, "<- Csc2updown_new \n");
#endif
}

/******************************************************************************/
/* void CscDiagDom(CscMatrix *cscmtx)                                         */
/*                                                                            */
/* Transforme la csc en csc a diagonale dominante                             */
/*                                                                            */
/* cscmtx : Csc                                                               */
/******************************************************************************/
/* !!! INUTILISEE, non testee !!! */
void CscDiagDom(CscMatrix *cscmtx)
{
  pastix_int_t itercblk;
  pastix_int_t itercoltab;
  pastix_int_t iterval;
  pastix_int_t itercol=0;
  pastix_int_t diag=0;
  double sum41;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscDiagDom \n");
#endif

  for (itercblk=0; itercblk<CSC_FNBR(cscmtx); itercblk++)
    {
      for (itercoltab=0;
	   itercoltab<CSC_COLNBR(cscmtx,itercblk);
	   itercoltab++)
	{
	  sum41 = 0;
	  for (iterval=CSC_COL(cscmtx,itercblk,itercoltab);
	       iterval<CSC_COL(cscmtx,itercblk,itercoltab+1);
	       iterval++)
	    {
	      if (CSC_ROW(cscmtx,iterval) != itercol)
		sum41 += ABS_FLOAT(CSC_VAL(cscmtx,iterval));
	      else
		diag = iterval;
	    }
	  CSC_VAL(cscmtx,diag) = sum41+1;
	  
	  itercol++;
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscDiagDom \n");
#endif
}

/* !!! INUTILISEE, non testee !!! */
pastix_int_t CscStrucSym(CscMatrix *cscmtx)
{
  /* csc_fnbr = 1 */
  pastix_int_t itercblk=0;
  pastix_int_t itercoltab;
  pastix_int_t iterval;
  pastix_int_t iterval2;
  pastix_int_t result=1;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscStrucSym \n");
#endif

  for (itercoltab=0;
       itercoltab<CSC_COLNBR(cscmtx,itercblk);
       itercoltab++)
    {
      for (iterval=CSC_COL(cscmtx,itercblk,itercoltab);
	   iterval<CSC_COL(cscmtx,itercblk,itercoltab+1);
	   iterval++)
	{
	  if (CSC_ROW(cscmtx,iterval) != itercoltab)
	    {
	      pastix_int_t rowidx=CSC_ROW(cscmtx, iterval);
	      pastix_int_t flag=0;
	      
	      for(iterval2=CSC_COL(cscmtx,itercblk,rowidx);
		  iterval2<CSC_COL(cscmtx,itercblk,rowidx+1);
		  iterval2++)
		{
		  if (CSC_ROW(cscmtx,iterval2) == itercoltab)
		    {
		      flag=1;
		      break;
		    }
		}
	      result = result && flag;
	    }
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscStrucSym \n");
#endif
  
  return result;
}
