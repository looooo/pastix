/*
 * File: out.h
 *
 * Define output format string for PaStiX.
 *
 * Authors:
 *   Mathieu Faverge - faverge@labri.fr
 *   Xavier  LACOSTE - lacoste@labri.fr
 *   Pierre  RAMET   - ramet@labri.fr
 */
#ifndef _OUT_H_
#define _OUT_H_

#define OUT_STARPU_TP         " StarPU : Thread policy : %s\n"
#define OUT_STARPU_STP        " StarPU : No thread policy, setting thread policy to : %s\n"
#define OUT_ENTETE_LINE1      " +--------------------------------------------------------------------+\n"
#define OUT_ENTETE_LINE2      " +              PaStiX : Parallel Sparse matriX package               +\n"
#define OUT_ENTETE_LINE3      " +--------------------------------------------------------------------+\n"
#define OUT_MATRIX_SIZE       "  Matrix size                                   %ld x %ld\n"
#define OUT_NNZ               "  Number of nonzeros in A                       %ld\n"
#define OUT_OPT_HEAD1         " +--------------------------------------------------------------------+\n"
#define OUT_OPT_HEAD2         " +  Options                                                           +\n"
#define OUT_OPT_HEAD3         " +--------------------------------------------------------------------+\n"
#define OUT_OPT_VERS          "        Version             :                   %s\n"
#define OUT_OPT_SMP           "        SMP_SOPALIN         :                   %s\n"
#define OUT_OPT_MPI           "        VERSION MPI         :                   %s\n"
#define OUT_OPT_DSCD          "        PASTIX_DYNSCHED     :                   %s\n"
#define OUT_OPT_STATS         "        STATS_SOPALIN       :                   %s\n"
#define OUT_OPT_NAPA          "        NAPA_SOPALIN        :                   %s\n"
#define OUT_OPT_IRECV         "        TEST_IRECV          :                   %s\n"
#define OUT_OPT_ISEND         "        TEST_ISEND          :                   %s\n"
#define OUT_OPT_THCOM         "        THREAD_COMM         :                   %s\n"
#define OUT_OPT_FUN           "        THREAD_FUNNELED     :                   %s\n"
#define OUT_OPT_TAG           "        TAG                 :                   %s\n"
#define OUT_OPT_OOC           "        OUT_OF_CORE         :                   %s\n"
#define OUT_OPT_DIST          "        PASTIX_DISTRIBUTED  :                   %s\n"
#define OUT_OPT_FORCE         "        FORCE_CONSO         :                   %s\n"
#define OUT_OPT_RFOB          "        RECV_FANIN_OR_BLOCK :                   %s\n"
#define OUT_OPT_METIS         "        METIS               :                   %s\n"
#define OUT_OPT_SCOTCH        "        WITH_SCOTCH         :                   %s\n"
#define OUT_OPT_INT           "        INTEGER TYPE        :                   %s\n"
#define OUT_OPT_FLOAT         "        FLOAT TYPE          :                   %s %s\n"
#define OUT_OPT_END           " +--------------------------------------------------------------------+\n"

#define OUT_STEP_ORDER        " Ordering :                                    \n"
#define OUT_SYMGRAPH          "   Symmetrizing graph                          \n"
#define OUT_NODIAG            "   Removing diag                               \n"
#define OUT_ORDERINIT         "   Initiating ordering                         \n"
#define OUT_ORDER_METHOD      "   Algorithm used: %s\n"

#define OUT_STEP_FAX          " Symbolic Factorization :                     \n"
#define OUT_FAX_METHOD        "   Algorithm used: %4s                       \n"
#define OUT_GLOBAL_NNZL       "   Number of nonzeroes in L structure      %ld\n"
#define OUT_GLOBAL_FILLIN     "   Fill-in                                 %lf\n"
#define OUT_GLOBAL_THFLOPCNT  "   Number of theoretical flop            %.5g %cflops\n"
#define OUT_GLOBAL_RLFLOPCNT  "   Number of performed flop              %.5g %cflops\n"

#define OUT_STEP_KASS         " Kass :                                       \n"
#define OUT_STEP_BLEND        " Analyse :                                    \n"
#define OUT_STEP_NUMFACT_LU   " Numerical Factorization (LU) :\n"
#ifdef TYPE_COMPLEX
#  define OUT_STEP_NUMFACT_LLT  " Numerical Factorization (LLh) :\n"
#else
#  define OUT_STEP_NUMFACT_LLT  " Numerical Factorization (LLt) :\n"
#endif
#define OUT_STEP_NUMFACT_LDLT " Numerical Factorization (LDLt) :\n"
#define OUT_STEP_NUMFACT_LDLH " Numerical Factorization (LDLh) :\n"
#define OUT_STEP_SOLVE        " Solve :                                      \n"
#define OUT_STEP_REFF         " Reffinement :                                \n"
#define TIME_COMPUTE_ORDERING "   Time to compute ordering                     %.3g s\n"
#define OUT_CLUSTNBR          "   Number of cluster                            %ld\n"
#define OUT_PROCNBR           "   Number of processor per cluster              %ld\n"
#define OUT_THRDNBR           "   Number of thread number per MPI process      %ld\n"
#define OUT_BLEND_CHKSMBMTX   "   Check the symbol matrix                      \n"
#define OUT_BLEND_CHKSOLVER   "   Check the solver structure                   \n"
#define OUT_BLEND_ELIMGRAPH   "   Building elimination graph                   \n"
#define OUT_BLEND_ELIMGRAPH2  "   Re-Building elimination graph                \n"
#define OUT_BLEND_COSTMATRIX  "   Building cost matrix                         \n"
#define OUT_BLEND_ELIMTREE    "   Building elimination tree                    \n"
#define OUT_BLEND_TASKGRAPH   "   Building task graph                          \n"
#define OUT_BLEND_NBTASK      "   Number of tasks                              %ld  \n"
#define OUT_BLEND_DISTPART    "   Distributing partition \n"
#define TIME_TO_ANALYSE       "   Time to analyze                              %.3g s\n"
#define NNZERO_WITH_FILLIN_TH "   Number of nonzeros in factorized matrix      %ld\n"
#define NNZERO_WITH_FILLIN    "%d : Number of nonzeros (local block structure) %ld\n"
#define SOLVMTX_WITHOUT_CO    "%d : SolverMatrix size (without coefficients)   %.3g %s\n"
#define OUT_FILLIN_TH         "   Fill-in                                      %lg\n"
#define NUMBER_OP_LU          "   Number of operations (LU)                    %g\n"
#define NUMBER_OP_LLT         "   Number of operations (LLt)                   %g\n"
#define TIME_FACT_PRED        "   Prediction Time to factorize (%s) %.3g s\n"
#define OUT_COEFSIZE          "   Maximum coeftab size (cefficients)           %.3g %s\n"
#define OUT_REDIS_CSC         "   Redistributing user CSC into PaStiX distribution\n"
#define OUT_REDIS_RHS         "   Redistributing user RHS into PaStiX distribution\n"
#define OUT_REDIS_SOL         "   Redistributing solution into Users' distribution\n"
#define OUT2_SOP_BINITG       "   --- Sopalin : Allocation de la structure globale ---\n"
#define OUT2_SOP_EINITG       "   --- Fin Sopalin Init                             ---\n"
#define OUT2_SOP_TABG         "   --- Initialisation des tableaux globaux          ---\n"
#define OUT2_SOP_BINITL       "   --- Sopalin : Local structure allocation         ---\n"
#define OUT2_SOP_NOTBIND      "   --- Sopalin : Threads are NOT binded             ---\n"
#define OUT2_SOP_BIND         "   --- Sopalin : Threads are binded                 ---\n"
#define OUT2_FUN_STATS        "     - %3ld : Envois %5ld - Receptions %5ld          -\n"
#define OUT2_SOP_BSOP         "   --- Sopalin Begin                                ---\n"
#define OUT2_SOP_ESOP         "   --- Sopalin End                                  ---\n"
#define OUT4_UPDO_TIME_INIT   " [%d][%d] Solve initialization time : %lg s \n"
#define OUT4_UPDO_COMM_TIME   " [%d][%d] Solve communication time : %lg s \n"
#define OUT4_FACT_COMM_TIME   " [%d][%d] Factorization communication time : %lg s \n"
#define OUT2_SOP_DOWN         "   --- Down Step                                    ---\n"
#define OUT2_SOP_DIAG         "   --- Diag Step                                    ---\n"
#define OUT2_SOP_UP           "   --- Up Step                                      ---\n"
#define GEN_RHS_1             "   Generate RHS for X=1\n"
#define GEN_RHS_I             "   Generate RHS for X=i\n"
#define GEN_SOL_0             "   Generate X0=0\n"
#define OOC_MEM_LIM_PERCENT   "   OOC memory limit                             %d%% of needed (%.3g %s)\n"
#define OOC_MEM_LIM           "   OOC memory limit                             %.3g %s\n"
#define OOC_IN_STEP           "   [%2d] IN %s :\n"
#define OOC_WRITTEN           "   [%2d]   written                               %.3g %s, allocated : %.3g %s\n"
#define OOC_READ              "   [%2d]   read                                  %.3g %s\n"
#define OOC_ALLOCATED         "   [%2d]   Allocated                             %.3g %s\n"
#define OOC_MAX_ALLOCATED     "   [%2d]   Maximum allocated                     %.3g %s\n"
#define OUT_ITERRAFF_GMRES    "   GMRES :\n"
#define OUT_ITERRAFF_PIVOT    "   Simple refinement :\n"
#define OUT_ITERRAFF_BICGSTAB  "   BICGSTAB :\n"
#define OUT_ITERRAFF_GRAD     "   Conjuguate gradient :\n"
#define OUT_ITERRAFF_ITER     "    - iteration %d :\n"
#define OUT_ITERRAFF_TTS      "         time to solve                          %.3g s\n"
#define OUT_ITERRAFF_TTT      "         total iteration time                   %.3g s\n"
#define OUT_ITERRAFF_ERR      "         error                                  %.5g  \n"
#define OUT_ITERRAFF_NORMA    "         ||A||                                  %.5g  \n"
#define OUT_ITERRAFF_NORMR    "         ||r||                                  %.5g  \n"
#define OUT_ITERRAFF_NORMB    "         ||b||                                  %.5g  \n"
#define OUT_ITERRAFF_BDIVR    "         ||r||/||b||                            %.5g  \n"
#define OUT_REDISCSCDTIME     "   Time to redistribute cscd                    %.3g s\n"
#define OUT_FILLCSCTIME       "   Time to fill internal csc                    %.3g s\n"
#define OUT_MAX_MEM_AF_SOP    "   Max memory used after factorization          %.3g %s\n"
#define OUT_MEM_USED_AF_SOP   "   Memory used after factorization              %.3g %s\n"
#define MAX_MEM_AF_CL         "   Max memory used after clean                  %.3g %s\n"
#define MEM_USED_AF_CL        "   Memory used after clean                      %.3g %s\n"
#define OUT_STATIC_PIVOTING   "   Static pivoting                              %ld\n"
#define OUT_INERTIA           "   Inertia                                      %ld\n"
#define OUT_INERTIA_PIVOT     "   Inertia (NB: with pivoting)                  %ld\n"
#define OUT_ESP_NBTASKS       "   Number of tasks added by esp                 %ld\n"
#define OUT_TIME_FACT         "   Time to factorize                            %.3g s\n"
#define OUT_FLOPS_FACT        "   FLOPS during factorization                   %.5g %s\n"
#define OUT_TIME_SOLV         "   Time to solve                                %.3g s\n"
#define OUT_RAFF_ITER_NORM    "   Refinement                                   %ld iterations, norm=%.3g\n"
#define OUT_PREC1             "   ||b-Ax||/||b||                               %.3g\n"
#define OUT_PREC2             "   max_i(|b-Ax|_i/(|b| + |A||x|)_i              %.3g\n"
#define OUT_TIME_RAFF         "   Time for refinement                          %.3g s\n"
#define OUT_END               " +--------------------------------------------------------------------+\n"

/*
 * Macro to decide who prints information
 */
#define pastix_print(mpirank, thrdrank, fmt, ...) do { if((mpirank == 0) && (thrdrank == 0)) fprintf(stdout, fmt, ##__VA_ARGS__); } while(0)


#define MEMORY_WRITE(mem) ( ((mem) < 1<<10) ?                           \
                            ( (double)(mem) ) :                         \
                            ( ( (mem) < 1<<20 ) ?                       \
                              ( (double)(mem)/(double)(1<<10) ) :       \
                              ( ((mem) < 1<<30 ) ?                      \
                                ( (double)(mem)/(double)(1<<20) ) :     \
                                ( (double)(mem)/(double)(1<<30) ))))
#define MEMORY_UNIT_WRITE(mem) (((mem) < 1<<10) ?                       \
                                "o" :                                   \
                                ( ( (mem) < 1<<20 ) ?                   \
                                  "Ko" :                                \
                                  ( ( (mem) < 1<<30 ) ?                 \
                                    "Mo" :                              \
                                    "Go" )))

#define PRINT_FLOPS(flops) ( ((flops) < 1<<10) ?                        \
                             ( (double)(flops) ) :                      \
                             ( ( (flops) < 1<<20 ) ?                    \
                               ( (double)(flops)/(double)(1<<10) ) :    \
                               ( ((flops) < 1<<30 ) ?                   \
                                 ( (double)(flops)/(double)(1<<20) ) :  \
                                 ( (double)(flops)/(double)(1<<30) ))))
#define PRINT_FLOPS_UNIT(flops) ( ((flops) < 1<<10) ?                   \
                                  ( "FLOPS" ) :                         \
                                  ( ( (flops) < 1<<20 ) ?               \
                                    ( "KFLOPS" ) :                      \
                                    ( ((flops) < 1<<30 ) ?              \
                                      ( "MFLOPS" ) :                    \
                                      ( "GFLOPS" ))))

static inline double
printflopsv( double flops )
{
    static double ratio = (double)(1<<10);
    int unit = 0;

    while ( (flops > ratio) && (unit < 9) ) {
        flops /= ratio;
        unit++;
    }
    return flops;
}

static inline char
printflopsu( double flops )
{
    static char units[9] = { ' ', 'k', 'm', 'g', 't', 'p', 'e', 'z', 'y' };
    static double ratio = (double)(1<<10);
    int unit = 0;

    while ( (flops > ratio) && (unit < 9) ) {
        flops /= ratio;
        unit++;
    }
    return units[unit];
}

#endif /* _OUT_H_ */
