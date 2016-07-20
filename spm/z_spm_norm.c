/**
 * @file z_spm_norm.c
 *
 *  PaStiX spm computational routines.
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2015-06-01
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "spm.h"
#include "z_spm.h"
#include "frobeniusupdate.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_internal
 *
 * z_spmFrobeniusNorm - Compute the Frobenius norm of the non distributed given
 * spm structure.
 *
 *  ||A|| = sqrt( sum( a_ij ^ 2 ) )
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm from which the norm need to be computed.
 *
 *******************************************************************************
 *
 * @return
 *           The computed norm
 *
 *******************************************************************************/
double
z_spmFrobeniusNorm( const pastix_spm_t *spm )
{
    pastix_int_t i, j, baseval;
    double *valptr = (double*)spm->values;
    double scale = 1.;
    double sumsq = 0.;

    if (spm->mtxtype == PastixGeneral) {
        for(i=0; i <spm->nnzexp; i++, valptr++) {
            frobenius_update( 1, &scale, &sumsq, valptr );

#if defined(PRECISION_z) || defined(PRECISION_c)
            valptr++;
            frobenius_update( 1, &scale, &sumsq, valptr );
#endif
        }
    }
    else {
        pastix_int_t *colptr = spm->colptr;
        pastix_int_t *rowptr = spm->rowptr;
        int nb, dofi, dofj, ii, jj;
        pastix_int_t *dofs=spm->dofs;
        baseval = spmFindBase( spm );

        switch( spm->fmttype ){
        case PastixCSC:
            for(i=0; i<spm->n; i++, colptr++)
            {
                dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
                for(j=colptr[0]; j<colptr[1]; j++, rowptr++)
                {
                    dofj = ( spm->dof > 0 ) ? spm->dof : dofs[(*rowptr-baseval)+1] - dofs[(*rowptr-baseval)];
                    for(ii=0 ; ii < dofi; ii++)
                    {
                        for(jj=0; jj < dofj ;jj++,valptr++)
                        {
                            nb = ( i == (*rowptr-baseval) ) ? 1 : 2;
                            frobenius_update( nb, &scale, &sumsq, valptr );

#if defined(PRECISION_z) || defined(PRECISION_c)
                            valptr++;
                            frobenius_update( nb, &scale, &sumsq, valptr );
#endif
                        }
                    }
                }
            }
            break;
        case PastixCSR:
            for(i=0; i<spm->n; i++, rowptr++) {
                for(j=rowptr[0]; j<rowptr[1]; j++, colptr++, valptr++) {
                    nb = ( i == (*colptr-baseval) ) ? 1 : 2;
                    frobenius_update( nb, &scale, &sumsq, valptr );

#if defined(PRECISION_z) || defined(PRECISION_c)
                    valptr++;
                    frobenius_update( nb, &scale, &sumsq, valptr );
#endif
                }
            }
            break;
        case PastixIJV:
            for(i=0; i <spm->nnz; i++, valptr++, colptr++, rowptr++) {
                nb = ( *rowptr == *colptr ) ? 1 : 2;
                frobenius_update( nb, &scale, &sumsq, valptr );

#if defined(PRECISION_z) || defined(PRECISION_c)
                valptr++;
                frobenius_update( nb, &scale, &sumsq, valptr );
#endif
            }
            break;
        default:
            fprintf(stderr, "z_spmFrobeniusNorm: Unknown Format\n");
        }
    }

    return scale * sqrt( sumsq );
}

/**
 *******************************************************************************
 *
 * @ingroup spm_internal
 *
 * z_spmMaxNorm - Compute the Max norm of the non distributed given spm
 * structure.
 *
 *  ||A|| = max( abs(a_ij) )
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm from which the norm need to be computed.
 *
 *******************************************************************************
 *
 * @return
 *           The computed norm
 *
 *******************************************************************************/
double
z_spmMaxNorm( const pastix_spm_t *spm )
{
    pastix_int_t i;
    pastix_complex64_t *valptr = (pastix_complex64_t*)spm->values;
    double tmp, norm = 0.;

    for(i=0; i < spm->nnzexp; i++, valptr++) {
        tmp = cabs( *valptr );
        norm = norm > tmp ? norm : tmp;
    }
    return norm;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_internal
 *
 * z_spmInfNorm - Compute the Infinity norm of the non distributed given spm
 * structure given by the maximum column sum.
 *
 *  ||A|| = max_i( sum_j(|a_ij|) )
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm from which the norm need to be computed.
 *
 *******************************************************************************
 *
 * @return
 *           The computed norm
 *
 *******************************************************************************/
double
z_spmInfNorm( const pastix_spm_t *spm )
{
    pastix_int_t col, row, i, j, k, ii, jj, baseval;
    pastix_int_t dofi, dofj;
    pastix_int_t *dofs=spm->dofs;

    pastix_complex64_t *valptr = (pastix_complex64_t*)spm->values;
    double norm = 0.;
    double *sumcol;

    MALLOC_INTERN( sumcol, spm->gNexp, double );
    memset( sumcol, 0, spm->gNexp * sizeof(double) );
    baseval = spmFindBase( spm );

    switch( spm->fmttype ){
    case PastixCSC:
        /* Original */
        /*
        for( col=0; col < spm->gN; col++ )
        {
            for( i=spm->colptr[col]-baseval; i<spm->colptr[col+1]-baseval; i++ )
            {
                sumcol[spm->rowptr[i]-baseval] += cabs( valptr[i] );
            }
        }
         */

        /* Add the symmetric/hermitian part */
        /*
        if ( (spm->mtxtype == PastixHermitian) ||
             (spm->mtxtype == PastixSymmetric) )
        {
            for( col=0; col < spm->gN; col++ )
            {
                for( i=spm->colptr[col]-baseval+1; i<spm->colptr[col+1]-baseval; i++ )
                {
                    sumcol[col] += cabs( valptr[i] );
                }
            }
        }
         */

        /* Dofs */
        for( i=0; i < spm->n; i++)
        {
            dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            for(k=spm->colptr[i]; k < spm->colptr[i+1]; k++)
            {
                j = spm->rowptr[k - baseval] - baseval;
                dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
                row  = ( spm->dof > 0 ) ? j        : dofs[j];
                for(ii=0; ii < dofi; ii++)
                {
                    for(jj=0; jj < dofj; jj++, valptr++)
                    {
                        {
                            sumcol[row + jj] += cabs( *valptr );
                        }
                    }
                }
            }
        }
        /* Add the symmetric/hermitian part */
        if ( (spm->mtxtype == PastixHermitian) ||
             (spm->mtxtype == PastixSymmetric) )
        {
            valptr = (pastix_complex64_t*)spm->values;
            for(i=0; i < spm->n; i++)
            {
                col  = ( spm->dof > 0 ) ? i        : dofs[i];
                dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
                for(k=spm->colptr[i]; k < spm->colptr[i+1]; k++)
                {
                    j = spm->rowptr[k - baseval] - baseval;
                    dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
                    if(i != j)
                    {
                        for(ii=0; ii < dofi; ii++)
                        {
                            for(jj=0; jj < dofj; jj++, valptr++)
                            {
                                sumcol[col + ii] += cabs( *valptr );
                            }
                        }
                    }
                    else
                    {
                        valptr += dofi * dofj;
                    }
                }
            }
        }
        break;

    case PastixCSR:
        for( row=0; row < spm->gN; row++ )
        {
            for( i=spm->rowptr[row]-baseval; i<spm->rowptr[row+1]-baseval; i++ )
            {
                sumcol[spm->colptr[i]-baseval] += cabs( valptr[i] );
            }
        }

        /* Add the symmetric/hermitian part */
        if ( (spm->mtxtype == PastixHermitian) ||
             (spm->mtxtype == PastixSymmetric) )
        {
            for( row=0; row < spm->gN; row++ )
            {
                for( i=spm->rowptr[row]-baseval+1; i<spm->rowptr[row+1]-baseval; i++ )
                {
                    sumcol[row] += cabs( valptr[i] );
                }
            }
        }
        break;

    case PastixIJV:
        for(i=0; i < spm->nnz; i++)
        {
            sumcol[spm->colptr[i]-baseval] += cabs( valptr[i] );
        }

        /* Add the symmetric/hermitian part */
        if ( (spm->mtxtype == PastixHermitian) ||
             (spm->mtxtype == PastixSymmetric) )
        {
            for(i=0; i < spm->nnz; i++)
            {
                if(spm->rowptr[i] != spm->colptr[i])
                    sumcol[spm->rowptr[i]-baseval] += cabs( valptr[i] );
            }
        }
        break;

    default:
        memFree_null( sumcol );
        return PASTIX_ERR_BADPARAMETER;
    }

    for( i=0; i<spm->gNexp; i++) //gNexp
    {
        if(norm < sumcol[i])
        {
            norm = sumcol[i];
        }
    }
    memFree_null( sumcol );

    return norm;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_internal
 *
 * z_spmOneNorm - Compute the Oneinity norm of the non distributed given spm
 * structure fiven by the maximum row sum
 *
 *  ||A|| = max_j( sum_i(|a_ij|) )
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm from which the norm need to be computed.
 *
 *******************************************************************************
 *
 * @return
 *           The computed norm
 *
 *******************************************************************************/
double
z_spmOneNorm( const pastix_spm_t *spm )
{
    pastix_int_t col, row, i, j, k, ii, jj, baseval;
    pastix_int_t dofi, dofj;

    pastix_complex64_t *valptr = (pastix_complex64_t*)spm->values;
    double norm = 0.;
    double *sumrow;
    pastix_int_t* dofs=spm->dofs;

    MALLOC_INTERN( sumrow, spm->gNexp, double );
    memset( sumrow, 0, spm->gNexp * sizeof(double) );
    baseval = spmFindBase( spm );

    switch( spm->fmttype ){
    case PastixCSC:
        /*
         for( col=0; col < spm->gN; col++ )
         {
         for( i=spm->colptr[col]-baseval; i<spm->colptr[col+1]-baseval; i++ )
         {
         sumrow[col] += cabs( valptr[i] );
         }
         }
         */

        for(i=0; i < spm->n; i++)
        {
            col  = ( spm->dof > 0 ) ? i        : dofs[i];
            dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            for(k=spm->colptr[i]; k < spm->colptr[i+1]; k++)
            {
                j = spm->rowptr[k - baseval] - baseval;
                dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
                for(ii=0; ii < dofi; ii++)
                {
                    for(jj=0; jj < dofj; jj++, valptr++)
                    {
                        {
                            sumrow[col + ii] += cabs( *valptr );
                        }
                    }
                }
            }
        }
        /* Add the symmetric/hermitian part */
        if ( (spm->mtxtype == PastixHermitian) ||
             (spm->mtxtype == PastixSymmetric) )
        {
            valptr = (pastix_complex64_t*)spm->values;
            for(i=0; i < spm->n; i++)
            {
                dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
                for(k=spm->colptr[i]; k < spm->colptr[i+1]; k++)
                {
                    j = spm->rowptr[k - baseval] - baseval;
                    row  = ( spm->dof > 0 ) ? j        : dofs[j];
                    dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
                    if(i != j)
                    {
                        for(ii=0; ii < dofi; ii++)
                        {
                            for(jj=0; jj < dofj; jj++, valptr++)
                            {
                                sumrow[row + jj] += cabs( *valptr );
                            }
                        }
                    }
                    else
                    {
                        valptr += dofi * dofj;
                    }
                }
            }
        }
        break;

    case PastixCSR:
        for( row=0; row < spm->gN; row++ )
        {
            for( i=spm->rowptr[row]-baseval; i<spm->rowptr[row+1]-baseval; i++ )
            {
                sumrow[row] += cabs( valptr[i] );
            }
        }

        /* Add the symmetric/hermitian part */
        if ( (spm->mtxtype == PastixHermitian) ||
             (spm->mtxtype == PastixSymmetric) )
        {
            for( row=0; row < spm->gN; row++ )
            {
                for( i=spm->rowptr[row]-baseval+1; i<spm->rowptr[row+1]-baseval; i++ )
                {
                    sumrow[spm->colptr[i]-baseval] += cabs( valptr[i] );
                }
            }
        }
        break;

    case PastixIJV:
        for(i=0; i < spm->nnz; i++)
        {
            sumrow[spm->rowptr[i]-baseval] += cabs( valptr[i] );
        }

        /* Add the symmetric/hermitian part */
        if ( (spm->mtxtype == PastixHermitian) ||
             (spm->mtxtype == PastixSymmetric) )
        {
            for(i=0; i < spm->nnz; i++)
            {
                if(spm->rowptr[i] != spm->colptr[i])
                    sumrow[spm->colptr[i]-baseval] += cabs( valptr[i] );
            }
        }
        break;

    default:
        memFree_null( sumrow );
        return PASTIX_ERR_BADPARAMETER;
    }

    for( i=0; i<spm->gNexp; i++)
    {
        if(norm < sumrow[i])
        {
            norm = sumrow[i];
        }
    }
    memFree_null( sumrow );

    return norm;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * z_spmNorm - Compute the norm of an spm matrix
 *
 *******************************************************************************
 *
 * @param[in] type
 *          = PastixMaxNorm: Max norm
 *          = PastixOneNorm: One norm
 *          = PastixInfNorm: Infinity norm
 *          = PastixFrobeniusNorm: Frobenius norm
 *
 * @param[in] spm
 *          The spm structure describing the matrix.
 *
 *******************************************************************************
 *
 * @return
 *      \retval The norm of the matrix spm
 *
 *******************************************************************************/
double
z_spmNorm( int ntype,
           const pastix_spm_t *spm)
{
    double norm = 0.;

    if(spm == NULL)
    {
        return -1.;
    }

    switch( ntype ) {
    case PastixMaxNorm:
        norm = z_spmMaxNorm( spm );
        break;

    case PastixInfNorm:
        norm = z_spmInfNorm( spm );
        break;

    case PastixOneNorm:
        norm = z_spmOneNorm( spm );
        break;

    case PastixFrobeniusNorm:
        norm = z_spmFrobeniusNorm( spm );
        break;

    default:
        fprintf(stderr, "z_spmNorm: invalid norm type\n");
        return -1.;
    }

    return norm;
}
