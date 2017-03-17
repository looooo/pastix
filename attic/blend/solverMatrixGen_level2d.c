    /****************************************/
    /** Compute the information for the    **/
    /** Forward and Backward triangular    **/
    /** Solution                           **/
    /****************************************/

    /* Pour l'instant uniquement si on est en 1d */
    if (ctrl->level2D == 0)
    {
    pastix_int_t *clust_mask       = NULL;
    pastix_int_t *clust_first_cblk = NULL;
    pastix_int_t *clust_highest    = NULL;
    pastix_int_t *uprecvcblk       = NULL;
    pastix_int_t  flag, cursor2;

        /** The initial symbol matrix is not expanded **/
        nodenbr = 0;
        for(i=0;i<symbmtx->cblknbr;i++)
            nodenbr += (symbmtx->cblktab[i].lcolnum-symbmtx->cblktab[i].fcolnum+1) * dof;
        solvmtx->updovct.gnodenbr= nodenbr;
        /*fprintf(stderr," GNODENBR %ld \n", (long)solvmtx->updovct.gnodenbr);*/

        /** Build the browtabs for each diagonal block **/
        MALLOC_INTERN(solvmtx->updovct.cblktab, solvmtx->cblknbr,UpDownCblk);
        cursor = 0;
        MALLOC_INTERN(clust_mask,       ctrl->clustnbr, pastix_int_t);
        MALLOC_INTERN(clust_first_cblk, ctrl->clustnbr, pastix_int_t);
        MALLOC_INTERN(clust_highest,    ctrl->clustnbr, pastix_int_t);

        solvmtx->updovct.downmsgnbr = 0;

        for(i=0;i<symbmtx->cblknbr;i++)
        {
            pastix_int_t brownbr;
            pastix_int_t fbrownum, lbrownum;

            /*if(cbprtab[i] == clustnum)*/
            if(simuctrl->bloktab[symbmtx->cblktab[i].bloknum].ownerclust == clustnum)
            {
                /*** Compute the list of clusters in the BROW (each cluster is list one time) ***/
                bzero(clust_mask, sizeof(pastix_int_t)*ctrl->clustnbr);
                bzero(clust_first_cblk, sizeof(pastix_int_t)*ctrl->clustnbr);
                for(j=0;j<ctrl->clustnbr;j++)
                    /*clust_highest[j] = - simuctrl->cblknbr;*/
                    clust_highest[j] = -1;


                brownbr = 0;
                fbrownum = symbmtx->cblktab[i].brownum;
                lbrownum = symbmtx->cblktab[i+1].brownum;
                assert( fbrownum != -1 );
                assert( lbrownum != -1 );
                for(j=fbrownum; j<lbrownum; j++)
                {
                    pastix_int_t cluster;
                    pastix_int_t cblk;
                    pastix_int_t bloknum = symbmtx->browtab[ j ];
                    cluster = simuctrl->bloktab[bloknum].ownerclust;
                    cblk = symbmtx->bloktab[bloknum].lcblknm;
#ifdef DEBUG_M
                    ASSERT( ctrl->candtab[cblk].treelevel   <= 0,   MOD_BLEND);
                    ASSERT( ctrl->candtab[cblk].distrib     == D1,  MOD_BLEND);
                    ASSERT( simuctrl->tasktab[cblk].cblknum == cblk,MOD_BLEND);
#endif

                    /*if( ctrl->candtab[cblk].treelevel >= clust_highest[cluster] )*/
                    if( simuctrl->tasktab[cblk].prionum >= clust_highest[cluster] )
                    {
                        /*clust_highest[cluster] = ctrl->candtab[cblk].treelevel;*/
                        clust_highest[cluster] = simuctrl->tasktab[cblk].prionum;
                        clust_first_cblk[cluster] = cblk;
                    }

                    if(clust_mask[cluster] == 0)
                    {
                        clust_mask[cluster] = 1;
                        brownbr++;
                    }
                }

                solvmtx->updovct.cblktab[cursor].browprocnbr = brownbr;
                if(solvmtx->updovct.cblktab[cursor].browprocnbr>0)
                {
                    MALLOC_INTERN(solvmtx->updovct.cblktab[cursor].browproctab,
                                  solvmtx->updovct.cblktab[cursor].browprocnbr,
                                  pastix_int_t);
                }
                else
                    solvmtx->updovct.cblktab[cursor].browproctab = NULL;

                if(clust_mask[clustnum] == 1)
                    solvmtx->updovct.cblktab[cursor].msgnbr = brownbr-1;
                else
                    solvmtx->updovct.cblktab[cursor].msgnbr = brownbr;
                solvmtx->updovct.downmsgnbr += solvmtx->updovct.cblktab[cursor].msgnbr;

                /*** Alloc the vector that will contain the global cblknum with the max priority for each processor in browproctab  ***/
                if(solvmtx->updovct.cblktab[cursor].browprocnbr>0)
                {
                    MALLOC_INTERN(solvmtx->updovct.cblktab[cursor].browcblktab,
                                  solvmtx->updovct.cblktab[cursor].browprocnbr,
                                  pastix_int_t);
                }
                else
                    solvmtx->updovct.cblktab[cursor].browcblktab = NULL;

                brownbr = 0;
                for(j=0;j<ctrl->clustnbr;j++)
                    if(clust_mask[j] == 1)
                    {
                        solvmtx->updovct.cblktab[cursor].browproctab[brownbr]   = j;
                        solvmtx->updovct.cblktab[cursor].browcblktab[brownbr++] = clust_first_cblk[j];
                    }

                solvmtx->updovct.cblktab[cursor].ctrbnbr = lbrownum - fbrownum;
                cursor++;
            }
        }

        /********************************************************************/
        /*** Find the list of local blocks in front of a diagonal blocks ****/
        /********************************************************************/
        cursor  = 0;
        cursor2 = 0;
        MALLOC_INTERN(solvmtx->updovct.gcblk2list, symbmtx->cblknbr, pastix_int_t);
        for(i=0;i<symbmtx->cblknbr;i++)
        {
            pastix_int_t fbrownum, lbrownum;

            fbrownum = symbmtx->cblktab[i].brownum;
            lbrownum = symbmtx->cblktab[i+1].brownum;

            solvmtx->updovct.gcblk2list[i] = -1;
            flag = 0;
            for(j=fbrownum; j<lbrownum; j++) {
                pastix_int_t bloknum = symbmtx->browtab[ j ];
                if( simuctrl->bloktab[ bloknum ].ownerclust == clustnum)
                {
                    if(flag == 0)
                    {
                        flag = 1;
                        solvmtx->updovct.gcblk2list[i] = cursor;
                        cursor++;
                    }
                    cursor2++;
                }
            }
        }
        solvmtx->updovct.gcblk2listnbr = symbmtx->cblknbr;

        MALLOC_INTERN(solvmtx->updovct.listptr, cursor+1, pastix_int_t);
        solvmtx->updovct.listptrnbr = cursor+1;
        MALLOC_INTERN(solvmtx->updovct.listblok, cursor2, pastix_int_t);
        MALLOC_INTERN(solvmtx->updovct.listcblk, cursor2, pastix_int_t);

        cursor  = 0;
        cursor2 = 0;
        for(i=0;i<symbmtx->cblknbr;i++)
        {
            flag = 0;
            for(j=0; j<ctrl->egraph->verttab[i].innbr;j++)
                if( simuctrl->bloktab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]].ownerclust == clustnum)
                {
                    if(flag == 0)
                    {
                        solvmtx->updovct.listptr[cursor2] = cursor;
                        cursor2++;
                        flag = 1;
                    }
#ifdef OOC
                    {
                        pastix_int_t tmp1,tmp2, tmp3;
                        pastix_int_t iter;
                        tmp1 = bloklocalnum[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j] ];
                        tmp2 = cblklocalnum[ctrl->egraph->ownetab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]]];
                        for (iter = solvmtx->updovct.listptr[cursor2-1]; iter < cursor; iter ++)
                            if (solvmtx->tasktab[tmp2].prionum <
                                solvmtx->tasktab[solvmtx->updovct.listcblk[iter]].prionum )
                            {
                                /* No problem with using solvmtx->updovct.listcblk[iter]
                                 * during first loop, cursor = 0, we don't use it,
                                 * and we set first.
                                 */
                                tmp3 = solvmtx->updovct.listcblk[iter];
                                solvmtx->updovct.listcblk[iter] = tmp2;
                                tmp2 = tmp3;

                                tmp3 = solvmtx->updovct.listblok[iter];
                                solvmtx->updovct.listblok[iter] = tmp1;
                                tmp1 = tmp3;
                            }
                        solvmtx->updovct.listblok[cursor] = tmp1;
                        solvmtx->updovct.listcblk[cursor] = tmp2;

                    }
#else
                    solvmtx->updovct.listblok[cursor] = bloklocalnum[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j] ];
                    solvmtx->updovct.listcblk[cursor] = cblklocalnum[ctrl->egraph->ownetab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]]];
#endif
                    cursor++;
                }

        }
        solvmtx->updovct.listptr[cursor2] = cursor;
        solvmtx->updovct.listnbr = cursor;

        solvmtx->updovct.loc2globnbr = solvmtx->cblknbr;
        MALLOC_INTERN(solvmtx->updovct.loc2glob, solvmtx->cblknbr, pastix_int_t);
        for(i=0;i<symbmtx->cblknbr;i++)
            if(cblklocalnum[i] >= 0)
                solvmtx->updovct.loc2glob[cblklocalnum[i]] = i;

        memFree_null(clust_mask);
        memFree_null(clust_first_cblk);
        memFree_null(clust_highest);

        /***** Fill lblk2gcblk ******/
        solvmtx->updovct.gcblknbr = symbmtx->cblknbr;
        MALLOC_INTERN(solvmtx->updovct.lblk2gcblk, symbmtx->bloknbr, pastix_int_t);
        for(i=0;i<symbmtx->bloknbr;i++)
            if(simuctrl->bloktab[i].ownerclust == clustnum)
                solvmtx->updovct.lblk2gcblk[bloklocalnum[i]] = symbmtx->bloktab[i].fcblknm;

        /* Compute the number of messages to receive during backward substitution */
        MALLOC_INTERN(uprecvcblk, symbmtx->cblknbr, pastix_int_t);
        for(i=0;i<symbmtx->cblknbr;i++)
            uprecvcblk[i] = 0;
        for (i=0; i<solvmtx->bublnbr; i++)
            for (j=0; j < solvmtx->ttsknbr[i]; j++)
            {
                SolverBlok * solvblok;
                cblknum = solvmtx->tasktab[solvmtx->ttsktab[i][j]].cblknum;
                for (solvblok =  solvmtx->cblktab[cblknum+1].fblokptr-1;
                     solvblok >= solvmtx->cblktab[cblknum].fblokptr+1; solvblok--)
                    /* if the contribution is not local */
                    if (solvmtx->bloktab[k].fcblknm <= 0)
                        uprecvcblk[solvmtx->updovct.lblk2gcblk[k]] = 1;
            }
        solvmtx->updovct.upmsgnbr = 0;
        for(i=0;i<symbmtx->cblknbr;i++)
            solvmtx->updovct.upmsgnbr += uprecvcblk[i];
        memFree_null(uprecvcblk);

        /*********************************/
        /*     Temporaire                */
        /*  Pour tester descente remonte */
        /*********************************/
        build_smx(&(solvmtx->updovct), symbmtx, simuctrl, ctrl);

    }
    /*********************** END TRIANGULAR INFO BUILDING ******************************************/

