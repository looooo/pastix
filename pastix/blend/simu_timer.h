/**
 *
 * @file simu_timer.h
 *
 *  PaStiX simulation routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#ifndef _SIMU_TIMER_H_
#define _SIMU_TIMER_H_

typedef struct SimuTimer_s {
  double s;
    /*  double ms;*/
} SimuTimer;

static inline int
timerComp(SimuTimer *t1, SimuTimer *t2)
{
    /* Return (t1 < t2) */
    if(t1->s < t2->s)
        return 1;
    else
        return 0;
}

static inline void
timerAdd(SimuTimer *timer, double t)
{
    timer->s += t;
}

static inline double
timerVal(const SimuTimer *t)
{
    return t->s;
}

static inline void
timerSet(SimuTimer *timer, double t)
{
    timer->s = t;
}

static inline void
timerSetMax(SimuTimer *timer, double t)
{
    if ( t > timer->s ) {
        timer->s = t;
    }
}

#endif /* _SIMU_TIMER_H_ */
