
/**
 *******************************************************************************
 *
 * @ingroup pastix_internal
 *
 * frobenius_update - Update the couple (scale, sumsq) with one element when
 * computing the Froebnius norm.
 *
 * The frobenius norm is equal to scale * sqrt( sumsq ), this method allows to
 * avoid overflow in the sum square computation.
 *
 *******************************************************************************
 *
 * @param[in,out] scale
 *           On entry, the former scale
 *           On exit, the update scale to take into account the value
 *
 * @param[in,out] sumsq
 *           On entry, the former sumsq
 *           On exit, the update sumsq to take into account the value
 *
 * @param[in] value
 *          The value to integrate into the couple (scale, sumsq)
 *
 *******************************************************************************/
static inline void
#if defined(PRECISION_d) || defined(PRECISION_z)
frobenius_update( int nb, double *scale, double *sumsq, double *value )
{
    double absval = fabs(*value);
    double ratio;
    if ( absval != 0. ){
        if ( (*scale) < absval ) {
            ratio = (*scale) / absval;
            *sumsq = (double)nb + (*sumsq) * ratio * ratio;
            *scale = absval;
        } else {
            ratio = absval / (*scale);
            *sumsq = (*sumsq) + (double)nb * ratio * ratio;
        }
    }
}
#elif defined(PRECISION_s) || defined(PRECISION_c)
frobenius_update( int nb, float *scale, float *sumsq, float *value )
{
    float absval = fabs(*value);
    float ratio;
    if ( absval != 0. ){
        if ( (*scale) < absval ) {
            ratio = (*scale) / absval;
            *sumsq = (float)nb + (*sumsq) * ratio * ratio;
            *scale = absval;
        } else {
            ratio = absval / (*scale);
            *sumsq = (*sumsq) + (float)nb * ratio * ratio;
        }
    }
}
#endif

