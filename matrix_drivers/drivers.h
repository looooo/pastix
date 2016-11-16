/**
 * @file drivers.h
 *
 *  $COPYRIGHTS$
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#ifndef _DRIVER_H_
#define _DRIVER_H_

#include "spm.h"

void pastix_ex_getoptions(int argc, char **argv,
                          pastix_int_t *iparam, double *dparam,
                          pastix_driver_t *driver, char **filename );


#endif /* _DRIVER_H_ */
