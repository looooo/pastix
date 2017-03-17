/************************************************************/
/**                                                        **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the blend control structure.        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     08 sep 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ Structure containing the structure passed through the blend primitives +*/
typedef struct BlendCtrl_ {
    pastix_int_t    count_ops ;      /*+ Print costs in term of number of elementary operations            +*/
    pastix_int_t    debug ;          /*+ Make additional checks after each step                            +*/
    pastix_int_t    timer;           /*+ Print execution times                                             +*/
    pastix_int_t    ooc;             /*+ Enable the out-of-core version of Pastix (Option unused for now)  +*/
    pastix_int_t    ricar;           /*+ Enable the ILU(k) dedicated steps                                 +*/
    pastix_int_t    leader;          /*+ Leader for sequential tasks                                       +*/

    /* Proportional Mapping */
    pastix_int_t    allcand;         /*+ All processors are candidate for each cblk                        +*/
    pastix_int_t    nocrossproc;     /*+ Forbid a processor to be candidate in two
                                      different branches shared with different partners                    +*/
    pastix_int_t    costlevel;       /*+ Enable/disable computation and use of subtree cost                +*/

    /* Spliting options */
    pastix_int_t    blcolmin ;       /*+ Minimun number of columns for a good use of BLAS primitives       +*/
    pastix_int_t    blcolmax;        /*+ Maximum number of columns for a good use of BLAS primitives       +*/
    pastix_int_t    abs;             /*+ Adaptative block size:
                                      - 0, all block are cut to blcolmin
                                      - >0, try to make (ncand*abs) cblk                                   +*/
    pastix_int_t    updatecandtab;   /*+ Update the candtab array after splitting, otherwise each new cblk
                                         has the same properties as the original one                       +*/

    /* 2D */
    pastix_int_t    autolevel;       /*+ Level to shift 1D to 2D is automaticaly computed                  +*/
    pastix_int_t    level2D;         /*+ number of levels to treat with a 2D distribution                  +*/
    pastix_int_t    ratiolimit;
    pastix_int_t    blblokmin ;      /*+ Minimum blocking size in 2D distribution                          +*/
    pastix_int_t    blblokmax;       /*+ Maximum blocking size in 2D distribution                          +*/

    /* Architecture */
    pastix_int_t    clustnum;        /*+ Id of current MPI process                                         +*/
    pastix_int_t    clustnbr;        /*+ Number of MPI processes                                           +*/
    pastix_int_t    total_nbcores;   /*+ Total number of physical cores used for the simulation            +*/
    pastix_int_t    total_nbthrds;   /*+ Total number of threads used for the simulation                   +*/
    pastix_int_t    local_nbcores;   /*+ Local number of physical cores used by the current MPI process    +*/
    pastix_int_t    local_nbthrds;   /*+ Local number of threads used by the current MPI process           +*/
    pastix_int_t    local_nbctxts;   /*+ Local number of contexts (used for dynamic scheduler and runtimes)+*/
    pastix_int_t   *clust2smp;       /*+ clust2smp[i] = SMP node on which i_th MPI
                                         process is running, if multiple MPI processes per node            +*/
    pastix_int_t   *core2clust;      /*+ core2clust[i] = cluster owning the core i                         +*/

    pastix_int_t   *iparm;           /*+ In/Out Integer parameters +*/
    double         *dparm;           /*+ In/Out Float parameters   +*/

    //BubbleTree        *btree;         /*+ arbre de bulles +*/
    EliminTree        *etree;         /*+ the elimination tree                   +*/
    CostMatrix        *costmtx;       /*+ the cost bounded to each cblk and blok +*/
    Cand              *candtab;       /*+ processor candidate tab                +*/
    Queue             *lheap;         /*+ Use to order leaves                    +*/
    ExtendVectorINT   *intvec;        /*+ vector of pastix_int_t used by several routines.
                                       The aim of this variable is to avoid
                                       repetedly memAlloc and memFree call      +*/
    ExtendVectorINT   *intvec2;       /*+ Another one                            +*/
    FILE              *tracefile;
} BlendCtrl;

int  blendCtrlInit ( BlendCtrl *ctrl,
                     pastix_int_t  clustnum,
                     pastix_int_t  clustnbr,
                     pastix_int_t *iparam,
                     double       *dparam );

void blendCtrlExit (BlendCtrl *);

void getCommunicationCosts( const BlendCtrl *ctrl,
                            pastix_int_t clustsrc,
                            pastix_int_t clustdst,
                            pastix_int_t sync_comm_nbr,
                            double *startup,
                            double *bandwidth);

