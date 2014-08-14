/** File: z_utils.c
 *
 *  @version 1.0.0
 *  @author Mathieu Faverge
 *  @author Pierre Ramet
 *  @author Xavier Lacoste
 *  @date 2011-11-11
 *  @precisions normal z -> c d s
 */
#ifndef Z_UTILS_H
#define Z_UTILS_H
#ifdef FORCE_NOMPI
#  define EXCHANGE_AX {}
#  define EXCHANGE_NORME {}
#  define IF_RANK_0 if (1)
#else
#  define EXCHANGE_AX(ax)                                               \
    {                                                                   \
        pastix_complex64_t * EAX_ax_rcv;                                \
        EAX_ax_rcv = malloc(globn*sizeof(pastix_complex64_t));          \
        MPI_Allreduce(ax, EAX_ax_rcv, globn, MPI_DOUBLE_COMPLEX, MPI_SUM, \
                      MPI_COMM_WORLD);                                  \
        free(ax);                                                       \
        ax = EAX_ax_rcv;                                                \
    }

#  define EXCHANGE_NORME(norme1, norme2)                                \
    {                                                                   \
        pastix_complex64_t EN_norme1_rcv, EN_norme2_rcv;                \
        MPI_Allreduce(&norme1, &EN_norme1_rcv, 1, MPI_DOUBLE_COMPLEX,   \
                      MPI_SUM, MPI_COMM_WORLD);                         \
        norme1 = EN_norme1_rcv;                                         \
        MPI_Allreduce(&norme2, &EN_norme2_rcv, 1, MPI_DOUBLE_COMPLEX,   \
                      MPI_SUM, MPI_COMM_WORLD);                         \
        norme2 = EN_norme2_rcv;                                         \
    }
#  define IF_RANK_0  if (mpid == 0)
#endif

#define PRINT_RHS_REAL(st, rh, nn, rk, verbose)         \
    {                                                   \
        if (verbose >= 5) {                             \
            int PRHS_ii;                                \
            fprintf(stdout,"%s (Proc %d) : ",st, rk);   \
            for (PRHS_ii= 0; PRHS_ii< nn; PRHS_ii++)    \
                fprintf(stdout,"%.3g ",rh[PRHS_ii]);    \
            fprintf(stdout,"\n");                       \
        }                                               \
    }
#define PRINT_RHS_TYPE_COMPLEX(st, rh, nn, rk, verbose)         \
    {                                                   \
        if (verbose >= 5) {                             \
            int PRHS_ii;                                \
            fprintf(stdout,"%s (Proc %d) : ",st, rk);   \
            for (PRHS_ii= 0; PRHS_ii< nn; PRHS_ii++)    \
                fprintf(stdout,"(%.3g %.3g) ",          \
                        creal(rh[PRHS_ii]),             \
                        cimag(rh[PRHS_ii]));            \
            fprintf(stdout,"\n");                       \
        }                                               \
    }

#ifdef TYPE_COMPLEX
#  define PRINT_RHS PRINT_RHS_TYPE_COMPLEX
#else
#  define PRINT_RHS PRINT_RHS_REAL
#endif

#define CONJ_REAL(x)  (x)
#define CONJ_TYPE_COMPLEX(x)  conjf(x)
#define CONJ_DTYPE_COMPLEX(x) conj(x)
#ifdef TYPE_COMPLEX
#  ifdef PREC_DOUBLE
#    define CONJ CONJ_DTYPE_COMPLEX
#  else
#    define CONJ CONJ_TYPE_COMPLEX
#  endif
#else
#  define CONJ CONJ_REAL
#endif

#define CHECK_SOL(sol, rhs, nn, rk)                                     \
    {                                                                   \
        int CS_ii,CS_jj;                                                \
        ax = malloc(nn*sizeof(pastix_complex64_t));                     \
        memset(ax, 0, nn*sizeof(pastix_complex64_t));                   \
        if (iparm[IPARM_TRANSPOSE_SOLVE] == API_NO)                     \
        {                                                               \
            for (CS_ii= 0; CS_ii < nn; CS_ii++)                         \
            {                                                           \
                for (CS_jj = colptr[CS_ii]-1;                           \
                     CS_jj < colptr[CS_ii+1] - 1;                       \
                     CS_jj++)                                           \
                {                                                       \
                    ax[rows[CS_jj]-1] += values[CS_jj]*sol[CS_ii];      \
                    if ((MTX_ISSYM(type) == 1) &&                       \
                        (CS_ii != (rows[CS_jj]-1)))                     \
                    {                                                   \
                        ax[CS_ii] += values[CS_jj]*sol[rows[CS_jj]-1];  \
                    }                                                   \
                    if ((MTX_ISHER(type) == 1) &&                       \
                        (CS_ii != (rows[CS_jj]-1)))                     \
                    {                                                   \
                        ax[CS_ii] = ax[CS_ii] +                         \
                            CONJ(values[CS_jj])*sol[rows[CS_jj]-1];     \
                    }                                                   \
                }                                                       \
            }                                                           \
        }                                                               \
        else                                                            \
        {                                                               \
            for (CS_ii= 0; CS_ii < nn; CS_ii++)                         \
            {                                                           \
                for (CS_jj = colptr[CS_ii]-1;                           \
                     CS_jj < colptr[CS_ii+1] - 1;                       \
                     CS_jj++)                                           \
                {                                                       \
                    ax[CS_ii] += values[CS_jj]*sol[rows[CS_jj]-1];      \
                    if ((MTX_ISSYM(type) == 1) &&                       \
                        (CS_ii != (rows[CS_jj]-1)))                     \
                    {                                                   \
                        ax[rows[CS_jj]-1] = ax[rows[CS_jj]-1] +         \
                            values[CS_jj]*sol[CS_ii];                   \
                    }                                                   \
                    if ((MTX_ISHER(type) == 1) &&                       \
                        (CS_ii != (rows[CS_jj]-1)))                     \
                    {                                                   \
                        ax[rows[CS_jj]-1] = ax[rows[CS_jj]-1] +         \
                            CONJ(values[CS_jj])*sol[rows[CS_jj]-1];     \
                    }                                                   \
                }                                                       \
            }                                                           \
        }                                                               \
        norme1= norme2 = 0;                                             \
        for (CS_ii= 0; CS_ii < nn; CS_ii++)                             \
        {                                                               \
            norme1 += (double)((ax[CS_ii] -                             \
                                rhs[CS_ii])*CONJ(ax[CS_ii] - rhs[CS_ii])); \
            norme2 += (double)(rhs[CS_ii] * CONJ(rhs[CS_ii]));          \
        }                                                               \
        if (rk == 0)                                                    \
            fprintf(stdout, "Precision : ||ax -b||/||b|| = %.20lg\n",   \
                    sqrt(norme1/norme2));                               \
        free(ax);                                                       \
    }

#define CHECK_DIST_SOL(colptr2, rows2, values2, rhs2, ncol2,            \
                       loc2glob2, globn, rhssaved_g)                    \
    {                                                                   \
                                                                        \
        pastix_int_t   * CDS_glob2loc;                                  \
        pastix_complex64_t * CDS_ax;                                    \
        pastix_complex64_t   CDS_norme1, CDS_norme2;                    \
        pastix_int_t     CDS_j, CDS_iter;                               \
        pastix_complex64_t * CDS_sol_g, *CDS_sol_g_recv;                \
        CDS_glob2loc = malloc(globn*sizeof(pastix_int_t));              \
        for (CDS_iter = 0; CDS_iter < globn; CDS_iter++)                \
            CDS_glob2loc[CDS_iter] = -1;                                \
        for (CDS_iter = 0; CDS_iter < ncol2; CDS_iter++)                \
            CDS_glob2loc[loc2glob2[CDS_iter]-1] = CDS_iter+1;           \
                                                                        \
        CDS_sol_g = malloc(globn*sizeof(pastix_complex64_t));           \
        memset(CDS_sol_g, 0, globn*sizeof(pastix_complex64_t));         \
        CDS_sol_g_recv = malloc(globn*sizeof(pastix_complex64_t));      \
        for (CDS_iter= 0; CDS_iter < ncol2; CDS_iter++)                 \
        {                                                               \
            CDS_sol_g[loc2glob2[CDS_iter]-1] = rhs2[CDS_iter];          \
        }                                                               \
        MPI_Allreduce(CDS_sol_g, CDS_sol_g_recv, globn,                 \
                      MPI_DOUBLE_COMPLEX, MPI_SUM,                      \
                      MPI_COMM_WORLD);                                  \
        free(CDS_sol_g);                                                \
        CDS_sol_g =CDS_sol_g_recv;                                      \
                                                                        \
        CDS_ax = malloc(globn*sizeof(pastix_complex64_t));              \
        memset(CDS_ax, 0, globn*sizeof(pastix_complex64_t));            \
        if (iparm[IPARM_TRANSPOSE_SOLVE] == API_NO)                     \
        {                                                               \
            for (CDS_iter= 0; CDS_iter < ncol2; CDS_iter++)             \
            {                                                           \
                for (CDS_j = colptr2[CDS_iter]-1;                       \
                     CDS_j < colptr2[CDS_iter+1] - 1; CDS_j++)          \
                {                                                       \
                    CDS_ax[rows2[CDS_j]-1] +=                           \
                        values2[CDS_j]*CDS_sol_g[loc2glob2[CDS_iter]-1]; \
                    if ((MTX_ISSYM(type) == 1) &&                       \
                        (loc2glob2[CDS_iter]-1 != (rows2[CDS_j]-1)))    \
                    {                                                   \
                        CDS_ax[loc2glob2[CDS_iter]-1] += values2[CDS_j]* \
                            CDS_sol_g[rows2[CDS_j]-1];                  \
                    }                                                   \
                    if ((MTX_ISHER(type) == 1) &&                       \
                        (loc2glob2[CDS_iter]-1 != (rows2[CDS_j]-1)))    \
                    {                                                   \
                        CDS_ax[loc2glob2[CDS_iter]-1] =                 \
                            CDS_ax[loc2glob2[CDS_iter]-1] +             \
                            CONJ(values2[CDS_j])*                       \
                            CDS_sol_g[rows2[CDS_j]-1];                  \
                    }                                                   \
                }                                                       \
            }                                                           \
        }                                                               \
        else                                                            \
        {                                                               \
            for (CDS_iter= 0; CDS_iter < ncol2; CDS_iter++)             \
            {                                                           \
                for (CDS_j = colptr2[CDS_iter]-1;                       \
                     CDS_j < colptr2[CDS_iter+1] - 1; CDS_j++)          \
                {                                                       \
                    CDS_ax[loc2glob2[CDS_iter]-1] +=                    \
                        values2[CDS_j]*CDS_sol_g[rows2[CDS_j]-1];       \
                    if ((MTX_ISSYM(type) == 1) &&                       \
                        (loc2glob2[CDS_iter]-1 != (rows2[CDS_j]-1)))    \
                    {                                                   \
                        CDS_ax[rows2[CDS_j]-1] += values2[CDS_j]*       \
                            CDS_sol_g[loc2glob2[CDS_iter]-1];           \
                    }                                                   \
                    if ((MTX_ISHER(type) == 1) &&                       \
                        (loc2glob2[CDS_iter]-1 != (rows2[CDS_j]-1)))    \
                    {                                                   \
                        CDS_ax[rows2[CDS_j]-1] =                        \
                            CDS_ax[rows2[CDS_j]-1] +                    \
                            CONJ(values2[CDS_j])*                       \
                            CDS_sol_g[loc2glob2[CDS_iter]-1];           \
                    }                                                   \
                }                                                       \
            }                                                           \
        }                                                               \
        free(CDS_sol_g);                                                \
        EXCHANGE_AX(CDS_ax);                                            \
        CDS_norme1= CDS_norme2 = 0;                                     \
        for (CDS_iter= 0; CDS_iter < ncol2; CDS_iter++)                 \
        {                                                               \
            CDS_norme1 +=                                               \
                (double)( ( CDS_ax[loc2glob2[CDS_iter]-1] -             \
                            rhssaved_g[loc2glob2[CDS_iter]-1] )*        \
                          CONJ(CDS_ax[loc2glob2[CDS_iter]-1] -          \
                               rhssaved_g[loc2glob2[CDS_iter]-1]));     \
            CDS_norme2 +=                                               \
                (double)((rhssaved_g[loc2glob2[CDS_iter]-1])*           \
                         CONJ(rhssaved_g[loc2glob2[CDS_iter]-1]));      \
        }                                                               \
        EXCHANGE_NORME(CDS_norme1, CDS_norme2);                         \
        IF_RANK_0 {                                                     \
            fprintf(stdout, "Precision : ||ax -b||/||b|| = %.20lg\n",   \
                    sqrt(CDS_norme1/CDS_norme2));                       \
        }                                                               \
        free(CDS_ax);                                                   \
        free(CDS_glob2loc);                                             \
    }
#endif /* not Z_UTILS_H */
