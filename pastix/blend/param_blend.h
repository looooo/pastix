/************************************************************/
/**                                                        **/
/**   NAME       : param_blend.h                           **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the parameters of blend.            **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     08 sep 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ The parameters structure definition +*/

typedef struct BlendParam_ {
  char * hpf_filename;          /*+ file name for HPF distribution  +*/
  char * trace_filename;        /*+ file name for Paragraph traces  +*/
  char * ps_filename;           /*+ file name for matrix postscript +*/

  pastix_int_t    hpf;             /*+ gener an HPF distribution file                                    +*/
  pastix_int_t    tracegen ;       /*+ gener a simulated Paragraph execution trace file                  +*/
  pastix_int_t    ps ;             /*+ gener a post-script of symbol matrix and elimination tree         +*/
  pastix_int_t    assembly;        /*+ Gener the info structure needed to assemble                       +*/
  char * solvmtx_filename;/*+ name of solver matrix files (in sequential mode                   +*/
  pastix_int_t    count_ops ;      /*+ print costs in term of number of elementary operations            +*/
  pastix_int_t    debug ;          /*+ make some check at certains execution points                      +*/
  pastix_int_t    timer;           /*+ print execution time                                              +*/
  pastix_int_t    recover;         /*+ take acount of a recover time estimation for ftgt                 +*/
  pastix_int_t    blcolmin ;       /*+ minimun number of column for a good use of BLAS primitives        +*/
  pastix_int_t    blcolmax;
  pastix_int_t    blblokmin ;      /*+ size of blockage for a good use of BLAS primitives  in 2D distribution +*/
  pastix_int_t    blblokmax;
  pastix_int_t    abs;             /*+ adaptative block size: := 0 all block are cut to blcolmin else try to make (ncand*abs) column +*/
  pastix_int_t    leader;          /*+ Processor leader for not parallele task (ex: gener assembly1D     +*/
  pastix_int_t    allcand;         /*+ All processor are candidat for a splitted cblk                    +*/
  pastix_int_t    nocrossproc;     /*+ Crossing processor forbiden in the splitting phase                +*/
  pastix_int_t    split;           /*+ Split cblk during proportional mapping algorithm                  +*/
  pastix_int_t    forceE2;
  pastix_int_t    level2D;         /*+ number of level to treat with a 2D distribution                   +*/
  pastix_int_t    candcorrect;
  pastix_int_t    clusterprop;     /*+ Proportionnal mapping with clustering for upper layers            +*/
  pastix_int_t    costlevel;       /*+ Calcul du cout de chaque sous arbre dans candtab                  +*/
  pastix_int_t    autolevel;       /*+ Level to shift 1D to 2D is automaticly computed                   +*/
  pastix_int_t    smpnbr;          /*+ Number of smp node                                                +*/
  pastix_int_t    procnbr;         /*+ Number of physical processors in a smp node                       +*/
  double ratiolimit;
  pastix_int_t    dense_endblock;   /*+ Treat the square right lower part of the matrix as a dense matrix+*/
  pastix_int_t    ooc;              /*+ To use the out of core version of Pastix                         +*/
  pastix_int_t    ricar;            /*+ If set to 1 then use blend dedicated to ricar                    +*/
  double oocmemlimit;      /*+ limit of physical memory for ooc                                 +*/
  pastix_int_t   *iparm;            /*+ In/Out Integer parameters +*/
  double *dparm;           /*+ In/Out Float parameters   +*/
  pastix_int_t     n;               /*+ Size of the matrix        +*/
} BlendParam;




pastix_int_t      blendParamInit(BlendParam *);
void     blendParamExit(BlendParam *);
