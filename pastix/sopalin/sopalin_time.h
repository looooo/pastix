/* 
   File: sopalin_time.h
   
   Few macros to manage clocks in sopalin.

   macros: 
     SOPALIN_CLOCK_INIT  - Initiate clock and start it.
     SOPALIN_CLOCK_STOP  - Stop clock.
     SOPALIN_CLOCK_GET   - Get value from clock.
     SOPALIN_CLOCK_TRACE - Get clock relatively to tracing start.
 */
#define SOPALIN_CLOCK_INIT  {clockInit(&(thread_data->sop_clk));clockStart(&(thread_data->sop_clk));}
#define SOPALIN_CLOCK_STOP  {clockStop(&(thread_data->sop_clk));}
#define SOPALIN_CLOCK_GET   clockVal(&(thread_data->sop_clk))
#define COMM_CLOCK_INIT  clockInit(&(thread_data->sop_clk_comm))
#define COMM_CLOCK_START clockStart(&(thread_data->sop_clk_comm))
#define COMM_CLOCK_STOP  clockStop(&(thread_data->sop_clk_comm))
#define COMM_CLOCK_GET   clockVal(&(thread_data->sop_clk_comm))
#define SOPALIN_CLOCK_TRACE (clockGet() - (sopalin_data->timestamp))
/* #define SOPALIN_CLOCK_TRACE ((clockGet() - (thread_data->sop_clk).time[0])) */
/* #define SOPALIN_CLOCK_TRACE (clockGet()) */
/* #define SOPALIN_CLOCK_TRACE (MPI_Wtime()) */

#ifdef COMPUTE_ALLOC 

static Clock alloc_clk;

#define ALLOC_CLOCK_INIT {clockInit(&alloc_clk);clockStart(&alloc_clk);}
#define ALLOC_CLOCK_STOP {clockStop(&alloc_clk);}
#define ALLOC_CLOCK_GET  clockVal(&alloc_clk)

#endif /* COMPUTE_ALLOC */
