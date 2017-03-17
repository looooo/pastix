#if defined(PASTIX_WITH_STARPU)
    /************************************************************************/
    /*  Fill the halo information                                           */
    /************************************************************************/
    if ( ctrl->iparm[IPARM_STARPU] == API_YES) {
        pastix_int_t halocblk=1;
        pastix_int_t bloknbr=0;
        SolverCblk * hcblk;
        SolverBlok * hblok;
        /* gcblk2halo[gcblk] == 0 : gcblk not local nor in halo
         *                   >  0 : local cblk number
         *                   <  0 : -halo cblk number
         */
        solvmtx->gcblknbr = symbmtx->cblknbr;
        MALLOC_INTERN(solvmtx->gcblk2halo, symbmtx->cblknbr, pastix_int_t);
        memset(solvmtx->gcblk2halo, 0, symbmtx->cblknbr*sizeof(pastix_int_t));
        for(i=0;i<symbmtx->cblknbr;i++) {
            if (cblklocalnum[i] >= 0) {
                /* If i is a local cblk */
                solvmtx->gcblk2halo[i] = cblklocalnum[i]+1;
                for( j=symbmtx->cblktab[i].bloknum;
                     j<symbmtx->cblktab[i+1].bloknum;
                     j++) {
                    pastix_int_t dst_cblk, dst_bloc;
                    dst_cblk = symbmtx->bloktab[j].cblknum;
                    if (solvmtx->gcblk2halo[dst_cblk] == 0 &&
                        symbmtx->cblktab[dst_cblk+1].bloknum -
                        symbmtx->cblktab[dst_cblk].bloknum > 0) {
                        dst_bloc = symbmtx->cblktab[dst_cblk].bloknum;
                        if ( simuctrl->bloktab[ dst_bloc ].ownerclust != clustnum) {
                            /* i updates remote cblk, add dst_cblk to halo */
                            solvmtx->gcblk2halo[dst_cblk] = -(halocblk++);
                            bloknbr += symbmtx->cblktab[dst_cblk+1].bloknum -
                                symbmtx->cblktab[dst_cblk].bloknum;
                        }
                    }
                }
            } else {
                /* If i is not a local cblk */
                if (pastix_starpu_with_fanin() == API_NO) {
                    if (solvmtx->gcblk2halo[i] == 0 &&
                        symbmtx->cblktab[i+1].bloknum -
                        symbmtx->cblktab[i].bloknum > 0) {
                        for(j=symbmtx->cblktab[i].bloknum;
                            j<symbmtx->cblktab[i+1].bloknum;
                            j++) {
                            pastix_int_t dst_cblk, dst_bloc;
                            dst_cblk = symbmtx->bloktab[j].cblknum;
                            dst_bloc = symbmtx->cblktab[dst_cblk].bloknum;
                            if ( simuctrl->bloktab[ dst_bloc ].ownerclust ==
                                 clustnum) {
                                /* i updates local cblk add i to halo*/
                                solvmtx->gcblk2halo[i] = -(halocblk++);
                                bloknbr += symbmtx->cblktab[i+1].bloknum -
                                    symbmtx->cblktab[i].bloknum;
                                break;
                            }
                        }
                    }
                }
            }
        }

        if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO) {
            fprintf(stdout, "%ld: Cblk number %ld, %ld blocks\n",
                        (long)solvmtx->clustnum, (long)solvmtx->cblknbr, (long)solvmtx->bloknbr);
            }
        if (pastix_starpu_with_fanin() == API_YES ) {
            pastix_int_t iter;
            pastix_int_t ftgtCblkIdx = 0;
            pastix_int_t ftgtBlokIdx;
            pastix_int_t clustnum;
            SolverCblk * fcblk;
            SolverBlok * fblok;
            MPI_Request * req;
            double fanin_coefnbr = 0;
            double fanin_coefnbr_pastix = 0;

            MALLOC_INTERN(solvmtx->fcblknbr, solvmtx->clustnbr, pastix_int_t);
            MALLOC_INTERN(solvmtx->fcblktab, solvmtx->clustnbr, SolverCblk*);
            MALLOC_INTERN(solvmtx->fbloktab, solvmtx->clustnbr, SolverBlok*);
            memset(solvmtx->fcblknbr, 0, solvmtx->clustnbr*sizeof(pastix_int_t));
            memset(solvmtx->fcblktab, 0, solvmtx->clustnbr*sizeof(SolverCblk*));
            memset(solvmtx->fbloktab, 0, solvmtx->clustnbr*sizeof(SolverBlok*));

            /**** OUTGOING FANIN ****/
            /* Count the number of Fanin blocks */
            for (ftgtBlokIdx = 0; ftgtBlokIdx < solvmtx->ftgtnbr; ftgtCblkIdx++) {
                FanInTarget * ftgt = &(solvmtx->ftgttab[ftgtBlokIdx]);
                pastix_int_t gcblk = ftgt->infotab[FTGT_GCBKDST];
                while( ftgtBlokIdx < solvmtx->ftgtnbr &&
                       ftgt->infotab[FTGT_GCBKDST] ==
                       gcblk) {
                    ftgtBlokIdx++;
                    ftgt++;
                }
            }

            if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO) {
                fprintf(stdout, "%ld: Outgoing Fanin cblk number %ld, %ld blocks\n",
                        (long)solvmtx->clustnum, (long)ftgtCblkIdx,
                        (long)solvmtx->ftgtnbr);

            }
            solvmtx->fcblknbr[solvmtx->clustnum]       = ftgtCblkIdx;
            MALLOC_INTERN(solvmtx->fcblktab[solvmtx->clustnum],
                          ftgtCblkIdx+1, SolverCblk);
            solvmtx->fbloktab[solvmtx->clustnum] = NULL;
            assert(ftgtBlokIdx == solvmtx->ftgtnbr);
            MALLOC_INTERN(solvmtx->fbloktab[solvmtx->clustnum],
                          solvmtx->ftgtnbr, SolverBlok);
            fcblk = solvmtx->fcblktab[solvmtx->clustnum];
            fblok = solvmtx->fbloktab[solvmtx->clustnum];
            /* Fill the outgoing fanin info */
            for (ftgtBlokIdx = 0; ftgtBlokIdx < solvmtx->ftgtnbr;) {
                FanInTarget * ftgt = &(solvmtx->ftgttab[ftgtBlokIdx]);
                fcblk->fcolnum = ftgt->infotab[FTGT_FCOLNUM];
                fcblk->lcolnum = ftgt->infotab[FTGT_LCOLNUM];
                fcblk->fblokptr = fblok;
                fcblk->stride  = 0;
                fcblk->procdiag = solvmtx->proc2clust[ftgt->infotab[FTGT_PROCDST]];
                fcblk->gcblknum = ftgt->infotab[FTGT_GCBKDST];
                /* While the target is the same we add bloks inside the
                 * fanin column block*/
                while( ftgtBlokIdx < solvmtx->ftgtnbr &&
                       ftgt->infotab[FTGT_GCBKDST] == fcblk->gcblknum) {
                    fblok->frownum = ftgt->infotab[FTGT_FROWNUM];
                    fblok->lrownum = ftgt->infotab[FTGT_LROWNUM];
                    fblok->coefind = fcblk->stride;
                    fcblk->stride +=
                        ftgt->infotab[FTGT_LROWNUM] -
                        ftgt->infotab[FTGT_FROWNUM] + 1;
                    fanin_coefnbr += (double)(cblk_colnbr(fcblk)*blok_rownbr(fblok));
                    fanin_coefnbr_pastix += (double)((ftgt->infotab[FTGT_LCOLNUM] -
                                                      ftgt->infotab[FTGT_FCOLNUM] + 1)
                                                     *blok_rownbr(fblok));

                    ftgtBlokIdx++;
                    fblok++;
                    ftgt++;
                }
                fcblk++;
            }

            if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO) {
                fprintf(stdout,
                        "%ld: Outgoing Fanin volume : %.3g coefficients (+%.3g%%),"
                        "                             %.3g with native scheduler\n",
                        (long)solvmtx->clustnum, fanin_coefnbr,
                        (fanin_coefnbr-fanin_coefnbr_pastix)/fanin_coefnbr_pastix,
                        fanin_coefnbr_pastix);
            }
            if (solvmtx->ftgtnbr > 0) {
                /*  virtual cblk to avoid side effect in the loops on cblk bloks */
                fcblk->fcolnum = (fcblk-1)->lcolnum+1;
                fcblk->lcolnum = (fcblk-1)->lcolnum+1;
                fcblk->fblokptr = fblok;
            } else {
                fcblk->fcolnum = 0;
                fcblk->lcolnum = 0;
                fcblk->fblokptr = solvmtx->fbloktab[solvmtx->clustnum];
            }
            assert(fcblk->fblokptr - solvmtx->fbloktab[solvmtx->clustnum] ==
                   solvmtx->ftgtnbr);

            /***** INCOMMING FANIN ***/
            for (clustnum = 0; clustnum<ctrl->clustnbr; clustnum++) {
                SimuCluster *simuclust = &(simuctrl->clustab[clustnum]);
                pastix_int_t ftgtnbr, fBlokNbr;
                if (clustnum == solvmtx->clustnum) continue;
                solvmtx->fcblknbr[clustnum] = 0;
                fBlokNbr = 0;

                /* Compute number of receiving contributions */
                ftgtnbr = extendint_Size(&(simuclust->ftgtsend[solvmtx->clustnum]));
                for(ftgtBlokIdx=0; ftgtBlokIdx<ftgtnbr;) {
                    pastix_int_t *infotab;
                    pastix_int_t gcblk;
                    ftgtnum = extendint_Read(&(simuclust->ftgtsend[solvmtx->clustnum]),
                                             ftgtBlokIdx);
                    infotab = simuctrl->ftgttab[ftgtnum].ftgt.infotab;
                    gcblk = infotab[FTGT_GCBKDST];
                    /* while still in same cblk go on */
                    while(ftgtBlokIdx<ftgtnbr &&
                          infotab[FTGT_GCBKDST] == gcblk) {
                        assert( solvmtx->proc2clust[infotab[FTGT_PROCDST]] ==
                                solvmtx->clustnum );

                        ftgtBlokIdx++;
                        fBlokNbr++;
                        ftgtnum =
                            extendint_Read(&(simuclust->ftgtsend[solvmtx->clustnum]),
                                           ftgtBlokIdx);
                        infotab = simuctrl->ftgttab[ftgtnum].ftgt.infotab;
                    }
                    solvmtx->fcblknbr[clustnum]++;
                }

                if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO) {
                    fprintf(stdout, "%ld: Fanin cblk number %ld,"
                            " %ld blocks received from %ld\n",
                            (long)solvmtx->clustnum, (long)solvmtx->fcblknbr[clustnum],
                            (long)fBlokNbr, (long)clustnum);

                }
                if(solvmtx->fcblknbr[clustnum] > 0) {
                    pastix_int_t ftgtnbr;
                    MALLOC_INTERN(solvmtx->fcblktab[clustnum],
                                  solvmtx->fcblknbr[clustnum]+1,
                                  SolverCblk);
                    MALLOC_INTERN(solvmtx->fbloktab[clustnum],
                                  fBlokNbr,
                                  SolverBlok);

                    fcblk = solvmtx->fcblktab[clustnum];
                    fblok = solvmtx->fbloktab[clustnum];

                    ftgtnbr = extendint_Size(&(simuclust->ftgtsend[solvmtx->clustnum]));
                    for(ftgtBlokIdx=0; ftgtBlokIdx<ftgtnbr;) {
                        pastix_int_t *infotab;
                        pastix_int_t ftgtnum;
                        ftgtnum = extendint_Read(&(simuclust->ftgtsend[solvmtx->clustnum]),
                                                 ftgtBlokIdx);
                        infotab = simuctrl->ftgttab[ftgtnum].ftgt.infotab;
                        fcblk->fcolnum = infotab[FTGT_FCOLNUM] * dof;
                        fcblk->lcolnum =
                            (infotab[FTGT_LCOLNUM] + 1) * dof - 1;
                        fcblk->fblokptr = fblok;
                        fcblk->stride  = 0;
                        fcblk->procdiag = solvmtx->clustnum;
                        fcblk->gcblknum = infotab[FTGT_GCBKDST];
                        while(ftgtBlokIdx < ftgtnbr &&
                              infotab[FTGT_GCBKDST] == fcblk->gcblknum) {
                            assert( solvmtx->proc2clust[infotab[FTGT_PROCDST]] ==
                                    solvmtx->clustnum );
                            assert( fblok - solvmtx->fbloktab[clustnum] < fBlokNbr );
                            fblok->frownum = infotab[FTGT_FROWNUM] * dof;
                            fblok->lrownum = (infotab[FTGT_LROWNUM]+1) *
                                dof - 1;
                            fblok->coefind = fcblk->stride;

                            fcblk->stride += infotab[FTGT_LROWNUM] -
                                infotab[FTGT_FROWNUM] + 1;
                            fblok++;
                            ftgtBlokIdx++;
                            ftgtnum =
                                extendint_Read(
                                    &(simuclust->ftgtsend[solvmtx->clustnum]),
                                    ftgtBlokIdx);
                            infotab = simuctrl->ftgttab[ftgtnum].ftgt.infotab;
                        }
                        fcblk++;
                    }

                    if (fcblk !=  solvmtx->fcblktab[clustnum]) {
                        /* virtual cblk to avoid side effect in the loops on
                         * cblk bloks */
                        fcblk->fcolnum = (fcblk-1)->lcolnum+1;
                        fcblk->lcolnum = (fcblk-1)->lcolnum+1;
                        fcblk->fblokptr = fblok;
                        fcblk->procdiag = -1;
                    } else {
                        fcblk->fcolnum = 0;
                        fcblk->lcolnum = 0;
                        fcblk->fblokptr = 0;
                        fcblk->procdiag = -1;
                    }
                }
            }
        }

        {
            /* Fill hcblktab and hbloktab */
            solvmtx->hcblknbr = halocblk-1;
            if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO) {
                fprintf(stdout, "%ld: Halo %ld cblks, %ld blocks\n",
                        (long)solvmtx->clustnum, (long)solvmtx->hcblknbr, (long)bloknbr);
            }

            MALLOC_INTERN(solvmtx->hcblktab, halocblk, SolverCblk);
            MALLOC_INTERN(solvmtx->hbloktab, bloknbr, SolverBlok);
            memset(solvmtx->gcblk2halo, 0, symbmtx->cblknbr*sizeof(pastix_int_t));

            hblok=solvmtx->hbloktab;
            hcblk=solvmtx->hcblktab;
            halocblk=0;
            for(i=0;i<symbmtx->cblknbr;i++) {
                if (cblklocalnum[i] >= 0) {
                    /* If i is a local cblk */
                    solvmtx->gcblk2halo[i] = cblklocalnum[i]+1;
                    for( j=symbmtx->cblktab[i].bloknum;
                         j<symbmtx->cblktab[i+1].bloknum;
                         j++) {
                        pastix_int_t dst_cblk, dst_bloc;
                        dst_cblk = symbmtx->bloktab[j].cblknum;
                        if (solvmtx->gcblk2halo[dst_cblk] == 0 &&
                            symbmtx->cblktab[dst_cblk+1].bloknum -
                            symbmtx->cblktab[dst_cblk].bloknum > 0) {
                            dst_bloc = symbmtx->cblktab[dst_cblk].bloknum;
                            if ( simuctrl->bloktab[ dst_bloc ].ownerclust !=
                                 clustnum) {
                                pastix_int_t bloc;
                                pastix_int_t coefind = 0;
                                /* i updates remote cblk, add dst_cblk to halo */
                                solvmtx->gcblk2halo[dst_cblk] = -(halocblk+1);
                                hcblk->fcolnum =
                                    symbmtx->cblktab[dst_cblk].fcolnum *
                                    dof;
                                hcblk->lcolnum =
                                    symbmtx->cblktab[dst_cblk].lcolnum *
                                    dof + dof-1;
                                hcblk->stride   = 0;
                                hcblk->fblokptr = hblok;
                                hcblk->procdiag = simuctrl->bloktab[
                                    dst_bloc ].ownerclust;
                                hcblk->gcblknum = dst_cblk;
                                for( bloc = symbmtx->cblktab[dst_cblk].bloknum;
                                     bloc < symbmtx->cblktab[dst_cblk+1].bloknum;
                                     bloc++) {
                                    pastix_int_t delta;
                                    delta  =  symbmtx->bloktab[bloc].lrownum -
                                        symbmtx->bloktab[bloc].frownum +1;
                                    delta *=  dof;
                                    hcblk->stride += delta;
                                    hblok->frownum = symbmtx->bloktab[bloc].frownum *
                                        dof;
                                    hblok->lrownum = symbmtx->bloktab[bloc].lrownum *
                                        dof + dof-1;
                                    hblok->coefind = coefind;
                                    coefind += delta;
                                    hblok ++;
                                }
                                halocblk++;
                                hcblk++;
                            }
                        }
                    }
                } else {
                    /* If i is not a local cblk */
                    if (pastix_starpu_with_fanin() == API_NO ) {

                        if (solvmtx->gcblk2halo[i] == 0 &&
                            symbmtx->cblktab[i+1].bloknum -
                            symbmtx->cblktab[i].bloknum > 0) {
                            for( j=symbmtx->cblktab[i].bloknum;
                                 j<symbmtx->cblktab[i+1].bloknum;
                                 j++) {
                                pastix_int_t dst_cblk, dst_bloc;
                                dst_cblk = symbmtx->bloktab[j].cblknum;
                                dst_bloc = symbmtx->cblktab[dst_cblk].bloknum;
                                if (simuctrl->bloktab[ dst_bloc ].ownerclust ==
                                    clustnum) {
                                    pastix_int_t bloc;
                                    pastix_int_t coefind = 0;
                                    /* i updates local cblk add i to halo*/
                                    solvmtx->gcblk2halo[i] = -(halocblk+1);
                                    hcblk->fcolnum =
                                        symbmtx->cblktab[i].fcolnum * dof;
                                    hcblk->lcolnum =
                                        symbmtx->cblktab[i].lcolnum * dof +
                                        dof-1;
                                    hcblk->stride   = 0;
                                    hcblk->fblokptr = hblok;
                                    hcblk->procdiag =
                                        simuctrl->bloktab[
                                            symbmtx->cblktab[i].bloknum ].ownerclust;
                                    hcblk->gcblknum = i;
                                    for( bloc = symbmtx->cblktab[i].bloknum;
                                         bloc < symbmtx->cblktab[i+1].bloknum;
                                         bloc++) {
                                        pastix_int_t delta;
                                        delta = symbmtx->bloktab[bloc].lrownum -
                                            symbmtx->bloktab[bloc].frownum +1;
                                        hcblk->stride += delta * dof;
                                        hblok->frownum =
                                            symbmtx->bloktab[bloc].frownum *
                                            dof;
                                        hblok->lrownum =
                                            symbmtx->bloktab[bloc].lrownum *
                                            dof + dof-1;
                                        hblok->coefind = coefind;
                                        coefind += delta;

                                        //hblok->cblknum = cblklocalnum[symbmtx->bloktab[j].cblknum];
                                        hblok++;
                                    }
                                    halocblk++;
                                    hcblk++;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            assert(halocblk == solvmtx->hcblknbr);
            if (halocblk > 0) {
                /*  virtual cblk to avoid side effect in the loops on cblk bloks */
                hcblk->fcolnum = solvmtx->hcblktab[halocblk-1].lcolnum+1;
                hcblk->lcolnum = solvmtx->hcblktab[halocblk-1].lcolnum+1;
                hcblk->fblokptr = hblok;
            } else {
                hcblk->fcolnum = 0;
                hcblk->lcolnum = 0;
                hcblk->fblokptr = hblok;
            }
            assert( bloknbr == hblok - solvmtx->hbloktab );
        }
        /************************************************************************/
        /*** Find the list of global column blocks contributing to each cblk ****/
        /************************************************************************/
        cursor  = 0;
        cursor2 = 0;
        MALLOC_INTERN(solvmtx->updovct.gcblk2glist, symbmtx->cblknbr, pastix_int_t);
        for(i=0;i<symbmtx->cblknbr;i++) {
            solvmtx->updovct.gcblk2glist[i] = -1;
            flag = 0;
            for(j=0; j<ctrl->egraph->verttab[i].innbr;j++) {
                if(flag == 0) {
                    flag = 1;
                    solvmtx->updovct.gcblk2glist[i] = cursor;
                    cursor++;
                }
                cursor2++;
            }
        }
        solvmtx->updovct.gcblk2glistnbr = symbmtx->cblknbr;

        MALLOC_INTERN(solvmtx->updovct.glistptr, cursor+1, pastix_int_t);
        solvmtx->updovct.glistptrnbr = cursor+1;
        MALLOC_INTERN(solvmtx->updovct.glistblok, cursor2, pastix_int_t);
        MALLOC_INTERN(solvmtx->updovct.glistcblk, cursor2, pastix_int_t);
        /* MALLOC_INTERN(solvmtx->updovct.glistproc, cursor2, pastix_int_t); */
        /* for (i = 0; i < cursor2; i++) solvmtx->updovct.glistproc[i] = -1; */

        cursor  = 0;
        cursor2 = 0;
        for(i=0;i<symbmtx->cblknbr;i++) {
            flag = 0;
            for(j=0; j<ctrl->egraph->verttab[i].innbr;j++) {
                pastix_int_t blocknum;
                if(flag == 0) {
                    solvmtx->updovct.glistptr[cursor2] = cursor;
                    cursor2++;
                    flag = 1;
                }
                blocknum = ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j];
                solvmtx->updovct.glistblok[cursor] = bloklocalnum[blocknum];
                solvmtx->updovct.glistcblk[cursor] = ctrl->egraph->ownetab[blocknum];
                /* solvmtx->updovct.glistproc[cursor] = proc2clust[ blprtab[blocknum] ]; */
                if (simuctrl->bloktab[blocknum].ownerclust == clustnum) {
                    ASSERT(solvmtx->gcblk2halo[ctrl->egraph->ownetab[blocknum]] > 0,
                           MOD_BLEND);
                } else {
                    ASSERT(solvmtx->gcblk2halo[ctrl->egraph->ownetab[blocknum]] <= 0,
                           MOD_BLEND);
                }
                cursor++;
            }
        }
        solvmtx->updovct.glistptr[cursor2] = cursor;
        solvmtx->updovct.glistnbr = cursor;
    }
#endif /* defined(PASTIX_WITH_STARPU) */
<<<<<<< a2d955e581fa3d7302669fbe075293a8e4c895d0

=======
>>>>>>> Further cleanup on blend
