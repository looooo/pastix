#ifndef SOPALIN_INIT_H
#define SOPALIN_INIT_H

/* Type of thread asking for allocation (must be a power of 2) */
#define INIT_COMPUTE   1
#define INIT_SEND      2
#define INIT_RECV      4
#define INIT_OOC       8

/************************************************/
/*     Structure pour le backup des données     */
/************************************************/

typedef struct Backup_ {
  pastix_int_t              cpftmax;         /* Double by LU version             */
  pastix_int_t              arftmax;         /* Double by LU version             */
  pastix_int_t              nbftmax;
  pastix_int_t *            task_ctrbcnt;    /* no inital value                  */
  pastix_int_t *            task_ftgtcnt;    /* no inital value                  */
  pastix_int_t *            fanin_ctrbnbr;   /* change updown information        */
  pastix_int_t *            fanin_prionum;   /* both used for tag and pack       */
  pastix_int_t *            btagtaskcnt;     /* btag taskcnt                     */
  pastix_int_t *            bcofsendcnt;     /* bcof sendcnt                     */
  pastix_int_t *            symbol_cblknum;  /* sopalin add negative information */
  pastix_int_t              symbol_nodenbr;  /* ???                              */
} Backup;

typedef struct BackupSolve_ {
  pastix_int_t *            fanin_ctrbnbr;   /* change updown information        */
  pastix_int_t *            symbol_cblknum;  /* sopalin add negative information */
} BackupSolve_t;


/* Allocate and initialize/Free globale data for solver */
void sopalin_init     (Sopalin_Data_t *sopalin_data, SolverMatrix *m, SopalinParam *sopaparam, int fact);
void sopalin_clean    (Sopalin_Data_t *sopalin_data, int step);

/* Allocate and initialize/Free thread data for solver */
void sopalin_init_smp (Sopalin_Data_t *sopalin_data, pastix_int_t me, int fact, int thrdtype);
void sopalin_clean_smp(Sopalin_Data_t *sopalin_data, pastix_int_t me);

/* Restore/backup des données modifiées pendant l'enchaînement facto/solve */
void sopalin_backup (SolverMatrix *datacode, Backup *b);
void sopalin_restore(SolverMatrix *datacode, Backup *b);

/* Restore/backup des données modifiées pendant le solve */
void solve_backup (SolverMatrix *datacode, BackupSolve_t *b);
void solve_restore(SolverMatrix *datacode, BackupSolve_t *b);

#endif /* SOPALIN_INIT_H */
