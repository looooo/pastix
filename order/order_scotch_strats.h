/**
 *
 * @file order_scotch_strats.h
 *
 *  PaStiX order routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains (PT-)Scotch strategy strings
 *
 * @version 5.1.0
 * @author Francois Pellegrini
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#ifndef _ORDER_SCOTCH_STRATS_H_
#define _ORDER_SCOTCH_STRATS_H_

#define NEW_STRATEGY
#if defined(NEW_STRATEGY)
#define SCOTCH_STRAT_DIRECT                                             \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?(m{rat=0.8,"                                 \
  ""                         "vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                      "m{rat=0.8,"                                  \
  ""                        "vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}};,"                         \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{rat=0.8,"                                 \
  ""                         "vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                        "m{rat=0.8,"                                \
  ""                          "vert=100,"                               \
  ""                          "low=h{pass=10},"                         \
  ""                          "asc=f{bal=0.2}};,"                       \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""      "ose=g}}"

  /* """unc=n{sep=/(levl<4)?(m{rat=0.8,"                                 \ */
  /* ""                         "vert=100,"                                \ */
  /* ""                         "low=h{pass=10},"                          \ */
  /* ""                         "asc=f{bal=0.2}})|"                        \ */
  /* ""                      "m{rat=0.8,"                                  \ */
  /* ""                        "vert=100,"                                 \ */
  /* ""                        "low=h{pass=10},"                           \ */
  /* ""                        "asc=f{bal=0.2}};,"                         \ */
  /* ""      "ole=f{cmin=10000,cmax=100000,frat=0.0},"                       \ */
  /* ""      "ose=g}}"                                                     \ */
#else
#define SCOTCH_STRAT_DIRECT                                             \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(levl<4)?(m{asc=j{move=200,pass=1000,bal=0.01},low=l{passhf=1,passgg=10,passch=1},type=h,vert=100,rat=0.8});," \
  ""      "ole=f{cmin=10000,cmax=100000,frat=0.0},"                       \
  ""      "ose=g},"    \
  """unc=n{sep=/(levl<4)?(m{asc=j{move=200,pass=1000,bal=0.01},low=l{passhf=1,passgg=10,passch=1},type=h,vert=100,rat=0.8});," \
  ""      "ole=f{cmin=10000,cmax=100000,frat=0.0},"                       \
  ""      "ose=g}}"

  /* "c{rat=0.7,"                                                          \ */
  /* """cpr=n{sep=/(levl<4)?(m{asc=j{move=200,pass=1000,bal=0.01},low=l{passhf=1,passgg=10,passch=1},type=h,vert=100,rat=0.8});," \ */
  /* ""      "ole=f{cmin=10000,cmax=100000,frat=0.0},"                       \ */
  /* ""      "ose=g},"    \ */
  /* """unc=n{sep=/(levl<4)?(m{asc=j{move=200,pass=1000,bal=0.01},low=l{passhf=1,passgg=10,passch=1},type=h,vert=100,rat=0.8});," \ */
  /* ""      "ole=f{cmin=10000,cmax=100000,frat=0.0},"                       \ */
  /* ""      "ose=g}}" */


  /* """unc=n{sep=/(levl<4)?(m{rat=0.8,"                                 \ */
  /* ""                         "vert=100,"                                \ */
  /* ""                         "low=l{passhf=1,passgg=1,passch=1},"                          \ */
  /* ""                         "asc=b{bnd=f{move=200,pass=1000,bal=0.01},org=(|l{passhf=1,passgg=10,passch=1})f{bal=0.2},width=1}})|"                        \ */
  /* ""                        "m{rat=0.8,"                                \ */
  /* ""                          "vert=100,"                               \ */
  /* ""                          "low=l{passhf=1,passgg=1,passch=1},"                         \ */
  /* ""                          "asc=b{bnd=f{move=200,pass=1000,bal=0.01},org=(|l{passhf=1,passgg=10,passch=1})f{bal=0.2},width=1}};,"                       \ */
  /* ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \ */
  /* ""      "ose=g}}" */

  /* """unc=n{sep=/(levl<4)?(m{asc=b{bnd=f{move=200,pass=1000,bal=0.01},org=(|l{passhf=1,passgg=10,passch=1})f{move=200,pass=1000,bal=0.01},width=3},low=l{passhf=1,passgg=10,passch=1},type=h,vert=100,rat=0.7});," \ */


  /* "n{sep=/(vert>120)?m{asc=b{bnd=f{move=200,pass=1000,bal=0.01},org=(|l{passhf=1,passgg=10,passch=1})f{move=200,pass=1000,bal=0.01},width=3},low=l{passhf=1,passgg=10,passch=1},type=h,vert=100,rat=0.7};,ole=s,ose=s}" */

/* "n{sep=/(levl<3)?m{asc=b{bnd=f{move=200,pass=1000,bal=0.01},org=(|l{passhf=1,passgg=10,passch=1})f{move=200,pass=1000,bal=0.01},width=3},low=l{passhf=1,passgg=10,passch=1},type=h,vert=100,rat=0.7};,ole=s,ose=s}" */
/* "n{sep=(/((vert)>(60))?((m{asc=b{bnd=f{move=200,pass=1000,bal=0.01},org=(|a{passdg=1,passgg=10})f{move=200,pass=1000,bal=0.01},width=3},low=a{passdg=1,passgg=10},type=h,vert=100,rat=0.7}));),ole=s,ose=s}" */
/* "n{sep=(/((vert)>(200))?((m{asc=b{bnd=j{move=200,pass=1000,bal=0.01},org=(|l{passhf=1,passgg=10,passch=1})j{move=200,pass=1000,bal=0.01},width=3},low=l{passhf=1,passgg=10,passch=1},type=h,vert=100,rat=0.7}));),ole=s,ose=s}" */
#endif

#define SCOTCH_STRAT_INCOMP                                             \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}}|"                          \
  ""                      "m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}};,"                         \
  ""      "ole=f{cmin=0,cmax=100000,frat=0.08},"                        \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                       "m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""            "ose=g}}"
#define SCOTCH_STRAT_PERSO                                              \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>%ld)?m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}}|"                          \
  ""                      "m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}};,"                         \
  ""      "ole=f{cmin=%ld,cmax=%ld,frat=%f},"                           \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>%ld)?(m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                       "m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=%ld,cmax=%ld,frat=%f},"                           \
  ""      "ose=g}}"

#define PTSCOTCH_STRAT_DIRECT                                           \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?m{rat=0.8,"                                 \
  ""                        "vert=100,"                                \
  ""                        "low=h{pass=10},"                          \
  ""                        "asc=f{bal=0.2}}|"                         \
  ""                      "m{rat=0.8,"                                 \
  ""                        "vert=100,"                                \
  ""                        "low=h{pass=10},"                          \
  ""                        "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=0,cmax=100000,frat=0.0},"                         \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{type=h,"                                  \
  ""                         "rat=0.8,"                                 \
  ""                         "vert=100000,"                             \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=08.2}})|"                       \
  ""                       "m{type=h,"                                  \
  ""                         "rat=0.8,"                                 \
  ""                         "vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""      "ose=g}}"
#define PTSCOTCH_STRAT_INCOMP                                           \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}}|"                          \
  ""                      "m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}};,"                         \
  ""      "ole=f{cmin=0,cmax=100000,frat=0.08},"                        \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                       "m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""      "ose=g}}"

#define PTSCOTCH_STRAT_PERSO  "c{rat=0.7,cpr=n{sep=/(vert>%ld)?m{vert=100,low=h{pass=10},asc=f{bal=0.2}}|m{vert=100,low=h{pass=10},asc=f{bal=0.2}};,ole=f{cmin=%ld,cmax=%ld,frat=%f},ose=g},unc=n{sep=/(vert>%ld)?(m{vert=100,low=h{pass=10},asc=f{bal=0.2}})|m{vert=100,low=h{pass=10},asc=f{bal=0.2}};,ole=f{cmin=%ld,cmax=%ld,frat=%f},ose=g}}"


#define SCOTCH_STRAT_CLIF                                               \
  "c{rat=0.7,"                                                          \
  "  cpr=n{sep=(/((levl)<(3))?((m{asc=b{bnd=f{move=200,pass=1000,bal=0.2},org=(|h{pass=10})f{move=200,pass=1000,bal=0.2},width=3},low=h{pass=10},type=h,vert=100,rat=0.7}|m{asc=b{bnd=f{move=200,pass=1000,bal=0.2},org=(|h{pass=10})f{move=200,pass=1000,bal=0.2},width=3},low=h{pass=10},type=h,vert=100,rat=0.7}));),ole=s,ose=s},unc=n{sep=(/((levl)<(3))?((m{asc=b{bnd=f{move=200,pass=1000,bal=0.2},org=(|h{pass=10})f{move=200,pass=1000,bal=0.2},width=3},low=h{pass=10},type=h,vert=100,rat=0.7}|m{asc=b{bnd=f{move=200,pass=1000,bal=0.2},org=(|h{pass=10})f{move=200,pass=1000,bal=0.2},width=3},low=h{pass=10},type=h,vert=100,rat=0.7}));),ole=s,ose=s}}"

#endif /* _ORDER_SCOTCH_STRATS_H_ */
