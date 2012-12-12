/**
 *
 * @file eztrace_convert_coreblas.c
 *
 *  PLASMA core_blas tracing kernels
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * This file provides the functions to generate the trace
 * in function of the events.
 *
 * @version 2.5.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 * Easy way to add new kernel:
 *   1 - Add the enum coreblas_ev_codes.h: COREBLAS_*kernelname* (Don't forget to check if it already exist or not)
 *   2 - Add the line to initialize your new kernel in the eztrace_convert_coreblas_init function
 *
 **/
#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <GTG.h>
#include <ev_codes.h>
#include <eztrace_list.h>
#include <eztrace_convert.h>
#include "coreblas_ev_codes.h"

#ifndef min
#define min( a, b ) ( (a) < (b) ? (a) : (b) )
#endif
#ifndef max
#define max( a, b ) ( (a) > (b) ? (a) : (b) )
#endif

#define COREBLAS_STATE      "ST_Thread"
#define COREBLAS_TASK_NAME  "Submitted Tasks counter"
#define COREBLAS_TASK_ALIAS "STasks"

#define COREBLAS_TASKR_NAME  "Global Ready Tasks counter"
#define COREBLAS_TASKR_ALIAS "GRTasks"

#define COREBLAS_TASKWR_NAME  "Local Ready Tasks counter"
#define COREBLAS_TASKWR_ALIAS "LRTasks"

#define COREBLAS_THREADS_MAX 4096


#define COREBLAS_NAME(ev) coreblas_array[ (int)( ((ev)->code) & COREBLAS_MASK_EVENTS) ].name
#define GTG_RANDOM gtg_get_random_color()

#define COREBLAS_INIT_EVENT( _idx_, _name_, _color_ )     \
    coreblas_array[_idx_].name  = _name_;                 \
    coreblas_array[_idx_].color = _color_;                \
    coreblas_array[_idx_].nb    = 0;                      \


/*
 * Check EZTrace version for compatibility
 */
#if !defined(EZTRACE_API_VERSION) || !(EZTRACE_API_VERSION > 0x00000400)
#error "EZTrace 0.7 or greater is required"
#endif

/*
 * Data for each event
 *   @name: name of the kernel (several kernels can use the same name, for
 *          example: gemm & gemm_f1, or laswp & laswpc
 *   @color: color assigne to the kernel
 *   @nb: number of calls to the kernel (with eztrace_stats)
 *   @sum: Total time spent executing this kernel (with eztrace_stats)
 *   @min: Minimal execution time of this kernel (with eztrace_stats)
 *   @max: Maximal execution time of this kernel (with eztrace_stats)
 */
typedef struct coreblas_s {
    char *name;
    gtg_color_t color;
    int  nb;
    double sum;
    double min;
    double max;
} coreblas_t;

static int         coreblas_array_initialized = 0;
static coreblas_t  coreblas_array[COREBLAS_NBMAX_EVENTS];
static gtg_color_t colors_array[20];

/*
 * Keep information on thread status.
 *   @tid: Thread Id
 *   @active: Thread is active/inactive
 *   @lasttime: Start time of this task
 */
typedef struct coreblas_thrdstate_s {
    unsigned int tid;
    int  active;
    double lasttime;
} coreblas_thrdstate_t;

static coreblas_thrdstate_t *thrdstate = NULL;
static int nbtrhd = 0;

static inline gtg_color_t gtg_get_random_color(){
    static int i = -1;
    i = (i+1)%20;
    return colors_array[i];
}

/*
 * Case with priority and request need to be handle correctly in GTG
 * before to be included here
 */
#ifdef TRACE_BY_SEQUENCE

#define MAX_SEQUENCE 100
typedef struct sequence_s {
    uint64_t id;
    char *name;
} sequence_t;

sequence_t *seqtab;

void sequenceInit(){
    seqtab = (sequence_t*)malloc(MAX_SEQUENCE * sizeof(sequence_t));
    memset(seqtab, 0, MAX_SEQUENCE * sizeof(sequence_t));
}

void sequenceDestroy(){
    int i=0;
    while( i < MAX_SEQUENCE && seqtab[i].id != 0)
    {
        free(seqtab[i].name);
        i++;
    }
    free(seqtab);
}

char *getSequence(uint64_t seq)
{
    int i=0;

    while ( (i < MAX_SEQUENCE)
            && (seqtab[i].id != 0)
            && (seqtab[i].id != seq) )
        i++;

    if (i < MAX_SEQUENCE)
    {
        if ( seqtab[i].id == seq )
        {
            return seqtab[i].name;
        }
        else
        {
            seqtab[i].id = seq;
            if ( asprintf(&(seqtab[i].name), "Sequence%03d", i) < 0 ) {
                fprintf(stderr, "Failed to create new sequence name\n");
                exit(-1);
            }

            addEntityValue(seqtab[i].name, COREBLAS_STATE, seqtab[i].name, colors_array[i%20] );
            return seqtab[i].name;
        }
    } else {
        fprintf(stderr, "WARNING: Too many sequences, you need to increase the limit and recompile\n");
        return "SequenceOutOfRange";
    }
}

void handle_coreblas_start (struct fxt_ev_64 *ev)
{
    FUNC_NAME;
    DECLARE_THREAD_ID_STR(_threadstr, CUR_INDEX, CUR_THREAD_ID);
    if ( GET_NBPARAMS(ev) > 0 ) {
        CHANGE() setState (CURRENT, COREBLAS_STATE, _threadstr, getSequence(GET_PARAM(ev, 2)) );
    } else {
        CHANGE() setState (CURRENT, COREBLAS_STATE, _threadstr, COREBLAS_NAME(ev) );
    }
}

#else

void handle_coreblas_start(struct fxt_ev_64 *ev)
{
    FUNC_NAME;
    DECLARE_THREAD_ID_STR(_threadstr, CUR_INDEX, CUR_THREAD_ID);
    if ( GET_NBPARAMS(ev) > 0 ) {
        CHANGE() setState (CURRENT, COREBLAS_STATE, _threadstr, COREBLAS_NAME(ev) );
    } else {
        CHANGE() setState (CURRENT, COREBLAS_STATE, _threadstr, COREBLAS_NAME(ev) );
    }
}

#endif

void handle_coreblas_task (struct fxt_ev_64 *ev)
{
    FUNC_NAME;
    DECLARE_PROCESS_ID_STR( process_id, CUR_INDEX );
    assert( GET_NBPARAMS(ev) == 1 );
    int value = (int)GET_PARAM(ev, 1);
    CHANGE() addVar (CURRENT, COREBLAS_TASK_ALIAS, process_id, (varPrec)value);
}

void handle_coreblas_taskw (struct fxt_ev_64 *ev)
{
    FUNC_NAME;
    assert( GET_NBPARAMS(ev) == 2 );
    DECLARE_PROCESS_ID_STR( process_id, CUR_INDEX );
    DECLARE_THREAD_ID_STR(  thread_id,  CUR_INDEX, (unsigned int)GET_PARAM(ev, 1));
    int value = (int)GET_PARAM(ev, 2);
    CHANGE() addVar (CURRENT, COREBLAS_TASKR_ALIAS,  process_id, (varPrec)value);
    CHANGE() addVar (CURRENT, COREBLAS_TASKWR_ALIAS, thread_id,  (varPrec)value);
}

void
handle_coreblas_stop ()
{
    FUNC_NAME;
    DECLARE_THREAD_ID_STR(_threadstr, CUR_INDEX, CUR_THREAD_ID);
    CHANGE() setState (CURRENT, COREBLAS_STATE, _threadstr, "wait");
}

int
eztrace_convert_coreblas_init()
{
    int i;

    if ( coreblas_array_initialized == 0 ) {

        /* Initialize the colors_array */
        colors_array[ 0] = GTG_RED;
        colors_array[ 1] = GTG_GREEN;
        colors_array[ 2] = GTG_BLUE;
        colors_array[ 3] = GTG_WHITE;
        colors_array[ 4] = GTG_TEAL;
        colors_array[ 5] = GTG_DARKGREY;
        colors_array[ 6] = GTG_YELLOW;
        colors_array[ 7] = GTG_PURPLE;
        colors_array[ 8] = GTG_LIGHTBROWN;
        colors_array[ 9] = GTG_DARKBLUE;
        colors_array[10] = GTG_PINK;
        colors_array[11] = GTG_DARKPINK;
        colors_array[12] = GTG_SEABLUE;
        colors_array[13] = GTG_KAKI;
        colors_array[14] = GTG_REDBLOOD;
        colors_array[15] = GTG_BROWN;
        colors_array[16] = GTG_GRENAT;
        colors_array[17] = GTG_ORANGE;
        colors_array[18] = GTG_MAUVE;
        colors_array[19] = GTG_LIGHTPINK;

        /* First initialization to fill in the gap */
        for(i=0; i<COREBLAS_NBMAX_EVENTS; i++) {
            coreblas_array[i].name  = "";
            coreblas_array[i].color = GTG_RANDOM;
            coreblas_array[i].nb    = -1;
            coreblas_array[i].sum   = 0.;
            coreblas_array[i].min   = 999999999999.;
            coreblas_array[i].max   = 0.;
        }

        /*
         * Some kernels have a predefined color to keep them from one figure to another.
         * Those which has never been assigned a color, use a random one
         */
        COREBLAS_INIT_EVENT(COREBLAS_GEMM,  "gemm",   GTG_YELLOW    );
        COREBLAS_INIT_EVENT(COREBLAS_HERK,  "herk",   GTG_WHITE     );
        COREBLAS_INIT_EVENT(COREBLAS_SYRK,  "syrk",   GTG_WHITE     );
        COREBLAS_INIT_EVENT(COREBLAS_HEMM,  "hemm",   GTG_DARKPINK  );
        COREBLAS_INIT_EVENT(COREBLAS_SYMM,  "symm",   GTG_DARKPINK  );
        COREBLAS_INIT_EVENT(COREBLAS_TRMM,  "trmm",   GTG_PURPLE    );
        COREBLAS_INIT_EVENT(COREBLAS_TRSM,  "trsm",   GTG_RED       );
        COREBLAS_INIT_EVENT(COREBLAS_HER2K, "her2k",  GTG_PINK      );
        COREBLAS_INIT_EVENT(COREBLAS_SYR2K, "syr2k",  GTG_PINK      );
        COREBLAS_INIT_EVENT(COREBLAS_GEMV,  "gemv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_GBMV,  "gbmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_HEMV,  "hemv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_HBMV,  "hbmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_HPMV,  "hpmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_SYMV,  "symv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_SBMV,  "sbmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_SPMV,  "spmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_TRMV,  "trmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_TBMV,  "tbmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_TPMV,  "tpmv",   GTG_TEAL      );
        COREBLAS_INIT_EVENT(COREBLAS_TRSV,  "trsv",   GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_TBSV,  "tbsv",   GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_TPSV,  "tpsv",   GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_GER,   "ger",    GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_GERU,  "geru",   GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_GERC,  "gerc",   GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_HER,   "her",    GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_HPR,   "hpr",    GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_HER2,  "her2",   GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_HPR2,  "hpr2",   GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_SYR,   "syr",    GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_SPR,   "spr",    GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_SYR2,  "syr2",   GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_SPR2,  "spr2",   GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_ROTG,  "rotg",   GTG_PURPLE    );
        COREBLAS_INIT_EVENT(COREBLAS_ROTMG, "rotmg",  GTG_PURPLE    );
        COREBLAS_INIT_EVENT(COREBLAS_ROT,   "rot",    GTG_PURPLE    );
        COREBLAS_INIT_EVENT(COREBLAS_ROTM,  "rotm",   GTG_PURPLE    );
        COREBLAS_INIT_EVENT(COREBLAS_SWAP,  "swap",   GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_SCAL,  "scal",   GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_COPY,  "copy",   GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_AXPY,  "axpy",   GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_GEADD, "geadd",  GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_DOT,   "dot",    GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_DOTU,  "dotu",   GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_DOTC,  "dotc",   GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_XDOT,  "xdot",   GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_NRM2,  "nrm2",   GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_ASUM,  "asum",   GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_AMAX,  "amax",   GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LACPY, "lacpy",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LANGE, "lange",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LANHE, "lanhe",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LANSY, "lansy",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LARFB, "larfb",  GTG_YELLOW    );
        COREBLAS_INIT_EVENT(COREBLAS_LARFT, "larft",  GTG_RED       );
        COREBLAS_INIT_EVENT(COREBLAS_LASWP, "laswp",  GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_LAUUM, "lauum",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_POTRF, "potrf",  GTG_GREEN     );
        COREBLAS_INIT_EVENT(COREBLAS_TRTRI, "trtri",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LASET, "laset",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LASSQ, "lassq",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_GELQT, "gelqt",  GTG_GREEN     );
        COREBLAS_INIT_EVENT(COREBLAS_GEQRT, "geqrt",  GTG_GREEN     );
        COREBLAS_INIT_EVENT(COREBLAS_GESSM, "gessm",  GTG_BLUE      );
        COREBLAS_INIT_EVENT(COREBLAS_GETRF, "getrf",  GTG_GREEN     );
        COREBLAS_INIT_EVENT(COREBLAS_LATRO, "latro",  GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_SSSSM, "ssssm",  GTG_YELLOW    );
        COREBLAS_INIT_EVENT(COREBLAS_TITRO, "titro",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_TRBMM, "trbmm",  GTG_BLUE      );
        COREBLAS_INIT_EVENT(COREBLAS_TRGMM, "trgmm",  GTG_BLUE      );
        COREBLAS_INIT_EVENT(COREBLAS_TSLQT, "tslqt",  GTG_RED       );
        COREBLAS_INIT_EVENT(COREBLAS_TSMLQ, "tsmlq",  GTG_YELLOW    );
        COREBLAS_INIT_EVENT(COREBLAS_TSMQR, "tsmqr",  GTG_YELLOW    );
        COREBLAS_INIT_EVENT(COREBLAS_TSQRT, "tsqrt",  GTG_RED       );
        COREBLAS_INIT_EVENT(COREBLAS_TSRFB, "tsrfb",  GTG_BLUE      );
        COREBLAS_INIT_EVENT(COREBLAS_TSTRF, "tstrf",  GTG_BLUE      );
        COREBLAS_INIT_EVENT(COREBLAS_TTLQT, "ttlqt",  GTG_REDBLOOD  );
        COREBLAS_INIT_EVENT(COREBLAS_TTMLQ, "ttmlq",  GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_TTMQR, "ttmqr",  GTG_ORANGE    );
        COREBLAS_INIT_EVENT(COREBLAS_TTQRT, "ttqrt",  GTG_REDBLOOD  );
        COREBLAS_INIT_EVENT(COREBLAS_TTRFB, "ttrfb",  GTG_SEABLUE   );
        COREBLAS_INIT_EVENT(COREBLAS_UNMLQ, "unmlq",  GTG_YELLOW    );
        COREBLAS_INIT_EVENT(COREBLAS_UNMQR, "unmqr",  GTG_YELLOW    );
        COREBLAS_INIT_EVENT(COREBLAS_GETRIP,"getrip", GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_PLGHE, "plghe",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_PLGSY, "plgsy",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_SHIFT, "shift",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_SHIFTW,"shiftw", GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_SWPAB, "swpab",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_PLRNT, "plrnt",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_PEMV,  "pemv",   GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_BRDALG,"brdalg", GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_TRDALG,"trdalg", GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_HEGST, "hegst",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_SYGST, "sygst",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_HERFB, "herfb",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_SYRFB, "syrfb",  GTG_RANDOM    );
        COREBLAS_INIT_EVENT(COREBLAS_LARFG,        "larfg",      GTG_BLUE     );
        COREBLAS_INIT_EVENT(COREBLAS_GEQP3_NORMS,  "qp3_norms",  GTG_YELLOW   );
        COREBLAS_INIT_EVENT(COREBLAS_GEQP3_PIVOT,  "qp3_pivot",  GTG_REDBLOOD );
        COREBLAS_INIT_EVENT(COREBLAS_GEQP3_UPDATE, "qp3_update", GTG_GREEN    );
        COREBLAS_INIT_EVENT(COREBLAS_SETVAR,       "setvar",     GTG_ORANGE   );

        coreblas_array_initialized = 1;
    }

    /*
     * Init only for trace conversion, not for statistics
     */
    if(get_mode() == EZTRACE_CONVERT) {
        addVarType( COREBLAS_TASK_ALIAS,   COREBLAS_TASK_NAME,   "CT_Process" );
        addVarType( COREBLAS_TASKR_ALIAS,  COREBLAS_TASKR_NAME,  "CT_Process" );
        addVarType( COREBLAS_TASKWR_ALIAS, COREBLAS_TASKWR_NAME, "CT_Thread"  );

        for(i=0; i<COREBLAS_NBMAX_EVENTS; i++) {
            if ( coreblas_array[i].nb == -1 )
                continue;

            addEntityValue( coreblas_array[i].name,
                            COREBLAS_STATE,
                            coreblas_array[i].name,
                            coreblas_array[i].color );
        }

        /* plasma coreblas */
        addEntityValue ("wait", COREBLAS_STATE, "wait", GTG_BLACK );

#ifdef TRACE_BY_SEQUENCE
        sequenceInit();
#endif
    }
    return 0;
}

int
handle_coreblas_events(struct fxt_ev_64 *ev)
{

    switch (ev->code)
        {
        case FUT_COREBLAS(STOP)  : handle_coreblas_stop();  break;
        case FUT_COREBLAS(TASK)  : handle_coreblas_task(ev);  break;
        case FUT_COREBLAS(TASKW) : handle_coreblas_taskw(ev); break;
        default:
            if ( ( (ev->code) & COREBLAS_PREFIX) ) {
                handle_coreblas_start(ev);
            }
            else
                return 0;
        }

    return 1;

}

void
eztrace_convert_coreblas_finalize()
{

}


int
handle_coreblas_stats(struct fxt_ev_64 *ev)
{
    int i;
    double time;

    if ( thrdstate == NULL ) {
        thrdstate = (coreblas_thrdstate_t*)malloc(COREBLAS_THREADS_MAX * sizeof(coreblas_thrdstate_t));
        memset( thrdstate, 0, COREBLAS_THREADS_MAX * sizeof(coreblas_thrdstate_t));
    }

    switch (ev->code) {
    case FUT_COREBLAS(STOP)  :
    {
        for (i=0; i<nbtrhd; i++) {
            if ( thrdstate[i].tid == (unsigned int)CUR_THREAD_ID) {
                if ( thrdstate[i].active == 0 ) {
                    fprintf(stderr, "WARNING: The end of a state appears before the beginning\n");
                    return 0;
                }

                time = ( CURRENT - thrdstate[i].lasttime );

                /* Check that we have an existing state */
                assert(  coreblas_array[ thrdstate[i].active ].nb >= 0 );

                if( coreblas_array[ thrdstate[i].active ].nb == 0 ) {
                    coreblas_array[ thrdstate[i].active ].sum = 0.;
                    coreblas_array[ thrdstate[i].active ].max = 0.;
                    coreblas_array[ thrdstate[i].active ].min = 999999999999.;
                }
                coreblas_array[ thrdstate[i].active ].nb++;
                coreblas_array[ thrdstate[i].active ].sum += time;
                coreblas_array[ thrdstate[i].active ].max = max( coreblas_array[ thrdstate[i].active ].max, time );
                coreblas_array[ thrdstate[i].active ].min = min( coreblas_array[ thrdstate[i].active ].min, time );

                thrdstate[i].active = 0;
                thrdstate[i].lasttime = 0;
                return 1;
            }
        }
        return 0;
    }
    break;

    case FUT_COREBLAS(TASK)  :
        break;
    case FUT_COREBLAS(TASKW) :
        break;

    default: /* All the different states */
        if ( ( (ev->code) & COREBLAS_PREFIX) ) {
            int ev_code = (int)((ev->code) & COREBLAS_MASK_EVENTS);
            for (i=0; i<nbtrhd; i++) {
                if ( thrdstate[i].tid == (unsigned int)CUR_THREAD_ID) {
                    if ( thrdstate[i].active != 0 ) {
                        fprintf(stderr, "WARNING: thread %d change to state %d before to stop previous state %d\n",
                                (int)CUR_THREAD_ID, ev_code, thrdstate[i].active );
                    }

                    thrdstate[i].active = ev_code;
                    thrdstate[i].lasttime = CURRENT;
                    return 1;
                }
            }

            /* Thread not found, we add it */
            if ( nbtrhd < COREBLAS_THREADS_MAX ) {
                thrdstate[nbtrhd].tid = (unsigned int)CUR_THREAD_ID;
                thrdstate[i].active = ev_code;
                thrdstate[nbtrhd].lasttime = CURRENT;
                nbtrhd++;
                return 1;
            }
        }
        return 0;
    }

    return 1;
}

/*
 * Print the results of statistics.
 */
void print_coreblas_stats() {
    int i;

    printf ( "\nCoreblas Module:\n");
    printf (   "-----------\n");

    for(i=0; i<COREBLAS_NBMAX_EVENTS; i++) {
        if ( coreblas_array[ i ].nb > 0 ) {
            printf ( "%s : %d calls\n"
                     "\tAverage time: %.3f ms\n"
                     "\tMaximun time: %.3f ms\n"
                     "\tMinimun time: %.3f ms\n",
                     coreblas_array[ i ].name,
                     coreblas_array[ i ].nb,
                     coreblas_array[ i ].sum / (double)(coreblas_array[ i ].nb),
                     coreblas_array[ i ].max,
                     coreblas_array[ i ].min);
        }
    }
}

struct eztrace_convert_module coreblas_module;

void libinit(void) __attribute__ ((constructor));
void libinit(void)
{
  coreblas_module.api_version = EZTRACE_API_VERSION;

  /* Specify the initialization function.
   * This function will be called once all the plugins are loaded
   * and the trace is started.
   * This function usually declared StateTypes, LinkTypes, etc.
   */
  coreblas_module.init = eztrace_convert_coreblas_init;

  /* Specify the function to call for handling an event
   */
  coreblas_module.handle = handle_coreblas_events;

  /* Specify the function to call for handling an event when eztrace_stats is called
   */
  coreblas_module.handle_stats = handle_coreblas_stats;

  /* Print the results of statistics
   */
  coreblas_module.print_stats = print_coreblas_stats;

  /* Specify the module prefix */
  coreblas_module.module_prefix = COREBLAS_EVENTS_ID;

  if ( asprintf(&coreblas_module.name, "coreblas") < 0 ) {
      fprintf(stderr, "Failed to create module name\n");
      exit(-1);
  }
  if ( asprintf(&coreblas_module.description, "Module for kernels used in PLASMA (BLAS, LAPACK and coreblas)") < 0 ) {
      fprintf(stderr, "Failed to create module description\n");
      exit(-1);
  }

  coreblas_module.token.data = &coreblas_module;

  /* Register the module to eztrace_convert */
  eztrace_convert_register_module(&coreblas_module);
}

void libfinalize(void) __attribute__ ((destructor));
void libfinalize(void)
{
#ifdef TRACE_BY_SEQUENCE
    if(get_mode() == EZTRACE_CONVERT) {
        sequenceDestroy();
    }
#endif
}
