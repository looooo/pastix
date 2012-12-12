/**
 *
 * @file plasmaos-hwloc.c
 *
 *  This file handles the mapping from pthreads calls to windows threads
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Piotr Luszczek
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef PLASMA_HWLOC

#include <hwloc.h>

static hwloc_topology_t plasma_topology = NULL; /* Topology object */
static volatile int     plasma_nbr = 0;

void plasma_topology_init(){

    pthread_mutex_lock(&mutextopo);
    if (!topo_initialized) {

        /* Allocate and initialize topology object.  */
        hwloc_topology_init(&plasma_topology);
        
        /* Perform the topology detection.  */
        hwloc_topology_load(plasma_topology);
        
        /* Get the number of cores (We don't want to use HyperThreading */
        sys_corenbr = hwloc_get_nbobjs_by_type(plasma_topology, HWLOC_OBJ_CORE);
        
        topo_initialized = 1;
    }
    plasma_nbr++;
    pthread_mutex_unlock(&mutextopo);
}

void plasma_topology_finalize(){

    plasma_unsetaffinity();
        
    pthread_mutex_lock(&mutextopo);
    plasma_nbr--;
    if ((topo_initialized ==1) && (plasma_nbr == 0)) {
        /* Destroy tpology */
        hwloc_topology_destroy(plasma_topology);
        
        topo_initialized = 0;
    }
    pthread_mutex_unlock(&mutextopo);
}

/**
 This routine will set affinity for the calling thread that has rank 'rank'.
 Ranks start with 0.

 If there are multiple instances of PLASMA then affinity will be wrong: all ranks 0
 will be pinned to core 0.

 Also, affinity is not restored when PLASMA_Finalize() is called, but is removed.
 */
int plasma_setaffinity(int rank) {
    hwloc_obj_t      obj;      /* Hwloc object    */ 
    hwloc_cpuset_t   cpuset;   /* HwLoc cpuset    */
    
    if (!topo_initialized) {
        plasma_error("plasma_setaffinity", "Topology not initialized");
        return PLASMA_ERR_UNEXPECTED;
    }

    /* Get last one.  */
    obj = hwloc_get_obj_by_type(plasma_topology, HWLOC_OBJ_CORE, rank);
    if (!obj)
        return PLASMA_ERR_UNEXPECTED;
    
    /* Get a copy of its cpuset that we may modify.  */
    /* Get only one logical processor (in case the core is SMT/hyperthreaded).  */
#if !defined(HWLOC_BITMAP_H)
    cpuset = hwloc_cpuset_dup(obj->cpuset);
    hwloc_cpuset_singlify(cpuset);
#else
    cpuset = hwloc_bitmap_dup(obj->cpuset);
    hwloc_bitmap_singlify(cpuset);
#endif
   
    /* And try to bind ourself there.  */
    if (hwloc_set_cpubind(plasma_topology, cpuset, HWLOC_CPUBIND_THREAD)) {
        char *str = NULL;
#if !defined(HWLOC_BITMAP_H)
        hwloc_cpuset_asprintf(&str, obj->cpuset);
#else
        hwloc_bitmap_asprintf(&str, obj->cpuset);
#endif
        printf("Couldn't bind to cpuset %s\n", str);
        free(str);
        return PLASMA_ERR_UNEXPECTED;
    }
    
    /* Get the number at Proc level ( We don't want to use HyperThreading ) */
    rank = obj->children[0]->os_index;
    
    /* Free our cpuset copy */
#if !defined(HWLOC_BITMAP_H)
    hwloc_cpuset_free(cpuset);
#else
    hwloc_bitmap_free(cpuset);
#endif
    return PLASMA_SUCCESS;
}

/**
 This routine will unset the affinity set by a previous call to 
 plasma_setaffinity.
 */
int plasma_unsetaffinity() {
    hwloc_obj_t      obj;      /* Hwloc object    */ 
    hwloc_cpuset_t   cpuset;   /* HwLoc cpuset    */
    
    if (!topo_initialized) {
        plasma_error("plasma_unsetaffinity", "Topology not initialized");
        return PLASMA_ERR_UNEXPECTED;
    }

    /* Get last one.  */
    obj = hwloc_get_obj_by_type(plasma_topology, HWLOC_OBJ_MACHINE, 0);
    if (!obj) {
        plasma_warning("plasma_unsetaffinity", "Could not get object");
        return PLASMA_ERR_UNEXPECTED;
    }
    
    /* Get a copy of its cpuset that we may modify.  */
    /* Get only one logical processor (in case the core is SMT/hyperthreaded).  */
#if !defined(HWLOC_BITMAP_H)
    cpuset = hwloc_cpuset_dup(obj->cpuset);
#else
    cpuset = hwloc_bitmap_dup(obj->cpuset);
#endif
   
    /* And try to bind ourself there.  */
    if (hwloc_set_cpubind(plasma_topology, cpuset, HWLOC_CPUBIND_THREAD)) {
        char *str = NULL;
#if !defined(HWLOC_BITMAP_H)
        hwloc_cpuset_asprintf(&str, obj->cpuset);
#else
        hwloc_bitmap_asprintf(&str, obj->cpuset);
#endif
        plasma_warning("plasma_unsetaffinity", "Could not bind to the whole machine");
        printf("Couldn't bind to cpuset %s\n", str);
        free(str);
        return PLASMA_ERR_UNEXPECTED;
    }
    
    /* Free our cpuset copy */
#if !defined(HWLOC_BITMAP_H)
    hwloc_cpuset_free(cpuset);
#else
    hwloc_bitmap_free(cpuset);
#endif
    return PLASMA_SUCCESS;
}

int plasma_getnuma_size() {
    hwloc_cpuset_t   cpuset;   /* HwLoc cpuset    */
    hwloc_obj_t      obj;
    int thrdnbr = 1;

    obj = hwloc_get_obj_by_type(plasma_topology, HWLOC_OBJ_NODE, 0);

    /* Get a copy of its cpuset that we may modify.  */
    if (obj != NULL) {
#if !defined(HWLOC_BITMAP_H)
      cpuset = hwloc_cpuset_dup(obj->cpuset);
#else
      cpuset = hwloc_bitmap_dup(obj->cpuset);
#endif
      thrdnbr = hwloc_get_nbobjs_inside_cpuset_by_type(plasma_topology, cpuset, HWLOC_OBJ_CORE);

      /* Free our cpuset copy */
#if !defined(HWLOC_BITMAP_H)
      hwloc_cpuset_free(cpuset);
#else
      hwloc_bitmap_free(cpuset);
#endif
    }

    return thrdnbr;
}
#ifdef __cplusplus
}
#endif

#endif /* PLASMA_HAS_COMPLEX */
