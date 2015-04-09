/**
 *
 * @file isched.h
 *
 * Copyright (c) 2008-2014 The University of Bordeaux, IPB, LaBRI, Inria -
 *                         Bordeaux-Sud-Ouest.  All rights reserved.
 *
 * Copyright (c) 2010-2014 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *  PaStiX thread binding routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains basic functions to bind threads.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#ifndef BINDTHREAD_H
#define BINDTHREAD_H

BEGIN_C_DECLS

enum isched_action_e {
    ISCHED_ACT_STAND_BY,
    ISCHED_ACT_PARALLEL,
    ISCHED_ACT_FINALIZE
};

#if defined(HAVE_HWLOC)
#include "isched_hwloc.h"
#define isched_topo_init               isched_hwloc_init
#define isched_topo_destroy            isched_hwloc_destroy
#define isched_topo_bind_on_core_index isched_hwloc_bind_on_core_index
#define isched_topo_unbind             isched_hwloc_unbind
#define isched_topo_world_size         isched_hwloc_world_size
#else
#define isched_topo_init               isched_nohwloc_init
#define isched_topo_destroy            isched_nohwloc_destroy
#define isched_topo_bind_on_core_index isched_nohwloc_bind_on_core_index
#define isched_topo_world_size         isched_nohwloc_world_size
#endif

int  isched_topo_init(void);
void isched_topo_destroy(void);
int  isched_topo_bind_on_core_index(int);
int  isched_topo_unbind();
int  isched_topo_world_size();

END_C_DECLS

#endif /* BINDTHREAD_H */
