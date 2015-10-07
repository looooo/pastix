/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#ifndef Z_UPDO_SENDRECV_H
#define Z_UPDO_SENDRECV_H

#define updo_up_WaitCtrb_storage API_CALL(updo_up_WaitCtrb_storage)
void updo_up_WaitCtrb_storage  ( z_Sopalin_Data_t *sopalin_data,
                                 pastix_int_t             updo_buffer_size,
                                 void           *updo_buffer,
                                 pastix_int_t             me,
                                 pastix_int_t             i);

#  define send_waitall      API_CALL(send_waitall)
#  define send_waitall_down API_CALL(send_waitall_down)
#  define send_waitall_up   API_CALL(send_waitall_up)
#  define send_waitone_down API_CALL(send_waitone_down)
#  define send_waitone_up   API_CALL(send_waitone_up)
#  define send_free_down    API_CALL(send_free_down)
#  define send_free_up      API_CALL(send_free_up)
#  define send_testall_down API_CALL(send_testall_down)
#  define send_testall_up   API_CALL(send_testall_up)

void send_free_down    ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t s_index );
void send_free_up      ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t s_index );
void send_testall_down ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me );
void send_testall_up   ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me );
void send_waitall      ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me,
                         void (*funcfree)(z_Sopalin_Data_t*, pastix_int_t, pastix_int_t) );
void send_waitall_down ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me );
void send_waitall_up   ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me );
int  send_waitone_down ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me );
int  send_waitone_up   ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me );

#ifdef PASTIX_UPDO_ISEND

#  define test_all_downsend API_CALL(test_all_downsend)
#  define test_all_upsend   API_CALL(test_all_upsend)
#  define wait_all_downsend API_CALL(wait_all_downsend)
#  define wait_all_upsend   API_CALL(wait_all_upsend)

void test_all_downsend ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me, int tag );
void test_all_upsend   ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me, int tag );
void wait_all_downsend ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me );
void wait_all_upsend   ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me );

#endif

/* Réceptions des contributions pour la remontée */
#define probe_updown               API_CALL(probe_updown)
#define updo_up_WaitCtrb_nostorage API_CALL(updo_up_WaitCtrb_nostorage)
#define updo_down_send             API_CALL(updo_down_send)
#define updo_down_recv             API_CALL(updo_down_recv)
#define updo_up_send               API_CALL(updo_up_send)
#define updo_up_recv               API_CALL(updo_up_recv)
#ifndef FORCE_NOMPI

/* /\* Thread de communication *\/ */
/* void* API_CALL(z_updo_thread_comm)(void * arg){ return NULL; } */
int  probe_updown ( MPI_Comm, pastix_int_t );

#  ifndef STORAGE
void updo_up_WaitCtrb_nostorage ( z_Sopalin_Data_t *sopalin_data,
                                  pastix_int_t, void *, pastix_int_t, pastix_int_t );
#  endif

/* Réceptions des comms MPI */
static inline
void updo_down_send ( z_Sopalin_Data_t *sopalin_data, pastix_int_t, pastix_int_t, pastix_int_t );
static inline
void updo_down_recv ( z_Sopalin_Data_t *sopalin_data, void *,
                      MPI_Status, pastix_int_t);
static inline
void updo_up_send   ( z_Sopalin_Data_t *sopalin_data, pastix_int_t, pastix_int_t, pastix_int_t);
static inline
void updo_up_recv   ( z_Sopalin_Data_t *sopalin_data, void *,
                      MPI_Status, pastix_int_t);

/* Thread de communication */
#  define z_updo_thread_comm API_CALL(z_updo_thread_comm)
void* z_updo_thread_comm ( void * );


#else /* FORCE_NOMPI */
void updo_up_WaitCtrb_nostorage ( z_Sopalin_Data_t *sopalin_data,
                                  pastix_int_t buf_size, void *buf,
                                  pastix_int_t me, pastix_int_t i)
{
  (void)sopalin_data; (void)buf_size; (void)buf; (void)me; (void)i;
}
int  probe_updown ( MPI_Comm comm, pastix_int_t tag )
{
  (void)comm; (void)tag;
  return 0;
}

#endif /* FORCE_NOMPI */

#endif /* not Z_UPDO_SENDRECV_H */
