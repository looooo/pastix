#define TRACE_PICL_TASK

#define TRACE_PICL_TASK_COMPUTE
/*
#define TRACE_PICL_TASK_COLOR
#define TRACE_PICL_TASK_COMPUTE
*/
/*
#define TRACE_PICL_TASK_SEND
#define TRACE_PICL_TASK_RECV
#define TRACE_PICL_COMM
*/
#ifdef TRACE_PICL_TASK_COMPUTE
#define PICL_DIAG     0
#define PICL_COMP_1D  1
#define PICL_E1       2
#define PICL_E2       3
#define PICL_ADD      4

#define PICL_SF       5
#define PICL_SB       6
#define PICL_RF       7
#define PICL_RB       8
#endif

#ifdef TRACE_PICL_TASK_COLOR
#define PICL_DIAG     (SOLV_COLOR(c))
#define PICL_COMP_1D  (SOLV_COLOR(c))
#define PICL_E1       (SOLV_COLOR(c))
#define PICL_E2       (SOLV_COLOR(c))
#define PICL_ADD      (SOLV_PROCNBR)

#define PICL_SF       (SOLV_PROCNBR+1)
#define PICL_SB       (SOLV_PROCNBR+2)
#define PICL_RF       (SOLV_PROCNBR+3)
#define PICL_RB       (SOLV_PROCNBR+4)
#endif

