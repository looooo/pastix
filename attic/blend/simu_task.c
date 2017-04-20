
double
simuTaskSendCost(SimuTask *taskptr, const pastix_int_t clustsrc, const pastix_int_t clustdst, BlendCtrl *ctrl)
{
    double startup, bandwidth;

    getCommunicationCosts( ctrl, clustsrc, clustdst,
                           ctrl->candtab[taskptr->cblknum].lccandnum -
                           ctrl->candtab[taskptr->cblknum].fccandnum + 1,
                           &startup, &bandwidth );

    assert( taskptr->taskid != COMP_1D );

    return (startup + bandwidth * taskptr->mesglen);
}
