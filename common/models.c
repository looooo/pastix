/**
 *
 * @file models.c
 *
 * PaStiX performance models routines
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @date 2017-11-09
 *
 **/
#include "common.h"
#include "models.h"

int
modelsGetKernelId( const char *kernelstr,
                   int        *nbcoef )
{
    if(0 == strcasecmp("getrf",  kernelstr)) { *nbcoef = 4; return PastixKernelGETRF; }
    if(0 == strcasecmp("hetrf",  kernelstr)) { *nbcoef = 4; return PastixKernelHETRF; }
    if(0 == strcasecmp("potrf",  kernelstr)) { *nbcoef = 4; return PastixKernelPOTRF; }
    if(0 == strcasecmp("pxtrf",  kernelstr)) { *nbcoef = 4; return PastixKernelPXTRF; }
    if(0 == strcasecmp("sytrf",  kernelstr)) { *nbcoef = 4; return PastixKernelSYTRF; }

    if(0 == strcasecmp("trsm1d", kernelstr)) { *nbcoef = 6; return PastixKernelTRSMCblk2d; }
    if(0 == strcasecmp("trsm2d", kernelstr)) { *nbcoef = 6; return PastixKernelTRSMBlok2d; }

    if(0 == strcasecmp("gemm1d", kernelstr)) { *nbcoef = 8; return PastixKernelGEMMCblk2d2d; }
    if(0 == strcasecmp("gemm2d", kernelstr)) { *nbcoef = 8; return PastixKernelGEMMBlok2d2d; }

    *nbcoef = 0;
    return -1;
}

pastix_model_t *
modelsNew()
{
    pastix_model_t *model = malloc(sizeof(pastix_model_t));

    int a, k;

    memset( model, 0, sizeof( pastix_model_t ) );

    for(a=0; a<4; a++) {
        for(k=0; k<PastixKernelLvl1Nbr; k++) {
            model->coefficients[a][k][0] = 0xdeadbeef;
        }
    }
    return model;
}

void
modelsFree( pastix_model_t *model )
{
    if ( model->name != NULL ) {
        free(model->name);
    }
    free(model->name);
}

void
modelsPropagate( pastix_model_t *model,
                 int arithm, pastix_ktype_t kernelid )
{
    double *coefs0 = model->coefficients[arithm][kernelid];
    double ratio;
    int a, k;
    int kstart = 0;
    int kend = -1;

    /* Look for loaded information about factorization kernels */
    if ( kernelid < PastixKernelSCALOCblk ) {
        for( k=PastixKernelGETRF; k<=PastixKernelSYTRF; k++) {
            if ( (k == (int)kernelid) || (model->coefficients[arithm][k][0] != 0xdeadbeef) ) {
                continue;
            }

            ratio = (( k == (int)PastixKernelGETRF ) ? 2. : 1. ) / (( kernelid == PastixKernelGETRF ) ? 2. : 1. );

            model->coefficients[arithm][k][0] =         coefs0[0];
            model->coefficients[arithm][k][1] =         coefs0[1];
            model->coefficients[arithm][k][2] = ratio * coefs0[2];
            model->coefficients[arithm][k][3] = ratio * coefs0[3];
        }

        for( a=0; a<4; a++) {
            if (a == arithm) {
                continue;
            }
            ratio = (0.5 * a + 0.5) / (0.5 * arithm + 0.5);

            for( k=PastixKernelGETRF; k<=PastixKernelSYTRF; k++) {
                if ( model->coefficients[a][k][0] != 0xdeadbeef ) {
                    continue;
                }

                model->coefficients[a][k][0] = ratio * model->coefficients[arithm][k][0];
                model->coefficients[a][k][1] = ratio * model->coefficients[arithm][k][1];
                model->coefficients[a][k][2] = ratio * model->coefficients[arithm][k][2];
                model->coefficients[a][k][3] = ratio * model->coefficients[arithm][k][3];
            }
        }
    }
    else if ( kernelid < PastixKernelTRSMCblk1d ) {
    }
    else if ( kernelid < PastixKernelGEMMCblk1d1d ) {
        kstart = PastixKernelTRSMCblk1d;
        kend   = PastixKernelTRSMBlok2d;
    }
    else {
        kstart = PastixKernelGEMMCblk1d1d;
        kend   = PastixKernelGEMMBlok2d2d;
    }

    for( k=kstart; k<=kend; k++) {
        if ( (k == (int)kernelid) || (model->coefficients[arithm][k][0] != 0xdeadbeef) ) {
            continue;
        }

        model->coefficients[arithm][k][0] = coefs0[0];
        model->coefficients[arithm][k][1] = coefs0[1];
        model->coefficients[arithm][k][2] = coefs0[2];
        model->coefficients[arithm][k][3] = coefs0[3];
        model->coefficients[arithm][k][4] = coefs0[4];
        model->coefficients[arithm][k][5] = coefs0[5];
        model->coefficients[arithm][k][6] = coefs0[6];
        model->coefficients[arithm][k][7] = coefs0[7];
    }

    for( a=0; a<4; a++) {
        if (a == arithm) {
            continue;
        }
        ratio = (0.5 * a + 0.5) / (0.5 * arithm + 0.5);

        for( k=kstart; k<=kend; k++) {
            if ( model->coefficients[a][k][0] != 0xdeadbeef ) {
                continue;
            }

            model->coefficients[a][k][0] = ratio * model->coefficients[arithm][k][0];
            model->coefficients[a][k][1] = ratio * model->coefficients[arithm][k][1];
            model->coefficients[a][k][2] = ratio * model->coefficients[arithm][k][2];
            model->coefficients[a][k][3] = ratio * model->coefficients[arithm][k][3];
            model->coefficients[a][k][4] = ratio * model->coefficients[arithm][k][4];
            model->coefficients[a][k][5] = ratio * model->coefficients[arithm][k][5];
            model->coefficients[a][k][6] = ratio * model->coefficients[arithm][k][6];
            model->coefficients[a][k][7] = ratio * model->coefficients[arithm][k][7];
        }
    }
}

int
modelsRead( pastix_model_t *model,
            const char *modelfilename )
{
    FILE *f = pastix_fopen( modelfilename );
    char *str, *strcoef;
    char kernelstr[7];
    int  rc, arithm, nbcoef;
    size_t strsize = 256;
    pastix_ktype_t kernelid;
    double *coefs;

    if ( f == NULL ) {
        fprintf(stderr, "Can't open model file\n");
        return -1;
    }

    str = malloc( strsize * sizeof(char) );
    do {
        rc = getline( &str, &strsize, f );
        if ( rc == -1 ) {
            perror( "modelsRead(getline header)" );
            return -1;
        }
    }
    while( str[0] == '#' );

    /* Read the model name */
    rc = getline( &str, &strsize, f );
    if ( rc == -1 ) {
        perror( "modelsRead(getline model name)" );
        return -1;
    }
    model->name = strdup( str );

    /* Read the model name */
    rc = getline( &str, &strsize, f );
    if ( rc == -1 ) {
        perror( "modelsRead(getline model name)" );
        return -1;
    }
    model->name = strdup( str );

    /* Read the model values */
    while( (rc = getline( &str, &strsize, f )) != 0 ) {
        if ( rc == -1 ) {
            perror( "modelsRead(getline values)" );
            return -1;
        }

        /* Skip commented lines */
        if ( str[0] == '#' ) {
            continue;
        }

        /* Read the arithmetic, and the kernel name */
        if ( sscanf( str, "%d;%6s;", &arithm, kernelstr ) != 2 ) {
            fprintf(stderr, "modelRead: Error reading line (%s)\n", str );
            continue;
        }

        if ( (arithm < 0) || (arithm > 3) ) {
            fprintf(stderr, "modelRead: Incorrect arithmetic in line (%s)\n", str );
            continue;
        }

        kernelid = modelsGetKernelId( kernelstr, &nbcoef );
        if ( (int)kernelid == -1 ) {
            fprintf(stderr, "modelRead: Incorrect kernel type in line (%s)\n", str );
            continue;
        }

        /* Read the corrrect number of coefficients and store them */
        coefs = model->coefficients[arithm][kernelid];
        strcoef = str + 3 + strlen( kernelstr );

        switch ( nbcoef ) {
        case 4:
            if ( sscanf( strcoef, "%le;%le;%le;%le",
                         coefs, coefs+1, coefs+2, coefs+3 ) != 4 )
            {
                fprintf(stderr, "modelRead: Pb reading the 4 coefficients in line (%s)\n", str );
                continue;
            }
            break;
        case 6:
            if ( sscanf( strcoef, "%le;%le;%le;%le;%le;%le",
                         coefs,   coefs+1, coefs+2,
                         coefs+3, coefs+4, coefs+5 ) != 6 )
            {
                fprintf(stderr, "modelRead: Pb reading the 6 coefficients in line (%s)\n", str );
                continue;
            }
            break;
        case 8:
            if ( sscanf( strcoef, "%le;%le;%le;%le;%le;%le;%le;%le",
                         coefs,   coefs+1, coefs+2, coefs+3,
                         coefs+4, coefs+5, coefs+6, coefs+7 ) != 8 )
            {
                fprintf(stderr, "modelRead: Pb reading the 8 coefficients in line (%s)\n", str );
                continue;
            }
            break;
        default:
            ;
        }

        modelsPropagate( model, arithm, kernelid );
    }

    fclose(f);
    free(str);

    return 0;
}

int
modelsInitDefaultCPU( pastix_model_t *model )
{
    int a = 1; /* Real double */
    int ktype;
    double *coefs;

    assert( model != NULL );

    /*
     * All coefficiensts given are for double real arithmetic
     */
    model->name = strdup("AMD1680 MKL");

    /* POTRF */
    ktype = PastixKernelPOTRF;
    coefs = model->coefficients[a][ktype];
    coefs[0] =  4.071507e-07;
    coefs[1] = -1.469893e-07;
    coefs[2] =  1.707006e-08;
    coefs[3] =  2.439599e-11;
    modelsPropagate( model, a, ktype );

    /* TRSM1D */
    ktype = PastixKernelTRSMCblk2d;
    coefs = model->coefficients[a][ktype];
    coefs[0] = 3.255168e-06;
    coefs[1] = 3.976198e-08;
    coefs[2] = 0.;
    coefs[3] = 0.;
    coefs[4] = 0.;
    coefs[5] = 2.626177e-10;
    modelsPropagate( model, a, ktype );

    /* GEMM1D */
    ktype = PastixKernelGEMMCblk2d2d;
    coefs = model->coefficients[a][ktype];
    coefs[0] =  1.216278e-06;
    coefs[1] =  0.;
    coefs[2] = -2.704179e-10;
    coefs[3] =  1.148989e-07;
    coefs[4] =  2.724804e-10;
    coefs[5] =  1.328900e-09;
    coefs[6] =  0.;
    coefs[7] =  2.429169e-10;
    modelsPropagate( model, a, ktype );

    return 0;
}

int
modelsInitDefaultGPU( pastix_model_t *model )
{
    int a = 1; /* Real double */
    int ktype;
    double *coefs;

    assert( model != NULL );

    /*
     * All coefficiensts given are for double real arithmetic
     */
    model->name = strdup("AMD1680 MKL");

    /* POTRF */
    ktype = PastixKernelPOTRF;
    coefs = model->coefficients[a][ktype];
    coefs[0] =  4.071507e-07;
    coefs[1] = -1.469893e-07;
    coefs[2] =  1.707006e-08;
    coefs[3] =  2.439599e-11;
    modelsPropagate( model, a, ktype );

    /* TRSM1D */
    ktype = PastixKernelTRSMCblk2d;
    coefs = model->coefficients[a][ktype];
    coefs[0] = 3.255168e-06;
    coefs[1] = 3.976198e-08;
    coefs[2] = 0.;
    coefs[3] = 0.;
    coefs[4] = 0.;
    coefs[5] = 2.626177e-10;
    modelsPropagate( model, a, ktype );

    /* GEMM1D */
    ktype = PastixKernelGEMMCblk2d2d;
    coefs = model->coefficients[a][ktype];
    coefs[0] =  1.216278e-06;
    coefs[1] =  0.;
    coefs[2] = -2.704179e-10;
    coefs[3] =  1.148989e-07;
    coefs[4] =  2.724804e-10;
    coefs[5] =  1.328900e-09;
    coefs[6] =  0.;
    coefs[7] =  2.429169e-10;
    modelsPropagate( model, a, ktype );

    return 0;
}

void
modelsLoad( pastix_data_t *pastix_data )
{
    char *filename = NULL;
    int rc = 0;

    /*
     * Get the model filename for the CPUs
     */
    pastix_data->cpu_models = modelsNew();
    filename = pastix_getenv( "PASTIX_MODELS_CPU" );

    if ( filename == NULL ) {
        rc = modelsInitDefaultCPU( pastix_data->cpu_models );
    }
    else {
        rc = modelsRead( pastix_data->cpu_models,
                         filename );
        pastix_cleanenv( filename );
    }
    if ( rc == -1 ) {
        modelsFree( pastix_data->cpu_models );
    }

    /*
     * Get the model filename for the GPUs
     */
    pastix_data->gpu_models = modelsNew();
    filename = pastix_getenv( "PASTIX_MODELS_GPU" );

    if ( filename == NULL ) {
        rc = modelsInitDefaultGPU( pastix_data->gpu_models );
    }
    else {
        rc = modelsRead( pastix_data->gpu_models,
                         filename );
        pastix_cleanenv( filename );
    }
    if ( rc == -1 ) {
        modelsFree( pastix_data->gpu_models );
    }
}
