#ifndef STARPU_SUBMIT_TASKS_H
#define STARPU_SUBMIT_TASKS_H

#define starpu_loop_data_  API_CALL(starpu_loop_data_)
#define starpu_loop_data_t API_CALL(starpu_loop_data_t)
typedef struct starpu_loop_data_ starpu_loop_data_t;

#define starpu_submit_one_trf API_CALL(starpu_submit_one_trf)
int starpu_submit_one_trf (pastix_int_t itertask, Sopalin_Data_t * sopalin_data);
#define starpu_submit_bunch_of_gemm API_CALL(starpu_submit_bunch_of_gemm)
int starpu_submit_bunch_of_gemm (pastix_int_t itertask, Sopalin_Data_t * sopalin_data);

#define starpu_submit_tasks API_CALL(starpu_submit_tasks)
int starpu_submit_tasks(Sopalin_Data_t   * sopalin_data);

#endif /* STARPU_SUBMIT_TASKS_H */
