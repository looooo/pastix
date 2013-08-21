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
    pastix_int_t    count_ops ;      /*+ Print costs in term of number of elementary operations            +*/
    pastix_int_t    debug ;          /*+ Make additional checks after each step                            +*/
    pastix_int_t    timer;           /*+ Print execution times                                             +*/
    pastix_int_t    ooc;             /*+ Enable the out-of-core version of Pastix (Option unused for now)  +*/
    pastix_int_t    ricar;           /*+ Enable the ILU(k) dedicated steps                                 +*/
    pastix_int_t    leader;          /*+ Leader for sequential tasks                                       +*/
    pastix_int_t    smpnbr;          /*+ Number of smp node                                                +*/
    pastix_int_t    procnbr;         /*+ Number of physical processors in a smp node                       +*/

    /* Proportional Mapping */
    pastix_int_t    allcand;         /*+ All processors are candidate for each cblk                        +*/
    pastix_int_t    nocrossproc;     /*+ Forbid a processor to be candidate in two
                                         different branches shared with different partners                 +*/
    pastix_int_t    costlevel;       /*+ Enable/disable computation and use of subtree cost                +*/

    /* Spliting options */
    pastix_int_t    blcolmin ;       /*+ Minimun number of columns for a good use of BLAS primitives       +*/
    pastix_int_t    blcolmax;        /*+ Maximum number of columns for a good use of BLAS primitives       +*/
    pastix_int_t    abs;             /*+ Adaptative block size:
                                           - 0, all block are cut to blcolmin
                                           - >0, try to make (ncand*abs) cblk                              +*/

    /* 2D */
    pastix_int_t    autolevel;       /*+ Level to shift 1D to 2D is automaticly computed                   +*/
    pastix_int_t    level2D;         /*+ number of level to treat with a 2D distribution                   +*/
    double          ratiolimit;
    pastix_int_t    blblokmin ;      /*+ Minimum blocking size in 2D distribution                          +*/
    pastix_int_t    blblokmax;       /*+ Maximum blocking size in 2D distribution                          +*/

    pastix_int_t   *iparm;           /*+ In/Out Integer parameters +*/
    double         *dparm;           /*+ In/Out Float parameters   +*/
} BlendParam;

pastix_int_t blendParamInit(BlendParam *, pastix_int_t, pastix_int_t *);
void         blendParamExit(BlendParam *);
