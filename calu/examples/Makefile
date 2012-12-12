###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.6
# @author Bilel Hadri
# @date 2010-11-15
#
###

PLASMA_DIR = ..
include ../Makefile.internal

ZSRC=	example_zgelqf.c  	\
	example_zgels.c 	\
	example_zgeqrs.c	\
	example_zgetrs.c  	\
	example_zpotrf.c  	\
	example_ztrsm.c   	\
        example_zgelqs.c  	\
	example_zgeqrf.c  	\
	example_zgesv.c   	\
	example_zposv.c   	\
	example_zpotrs.c  	\
	example_ztrsmpl.c	\
	example_zunmqr.c 	\
	example_zunmlq.c	\
	\
	example_zposv_f.f	\
	example_zgesv_f.f 	\
	example_zgels_f.f

-include Makefile.src

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ifeq (${PLASMA_F90}, 1)
ALLSRC:= $(ALLSRC)			\
	$(ZF90) $(CF90) $(DF90) $(SF90)	\
	example_zgetrf_tile_f90.f90	\
	example_cgeqrf_tile_f90.f90	\
	example_dgesv_tile_async_f90.f90\
	example_sgebrd_f90.f90
endif

ALLOBJ  = $(ALLSRC:.c=.o)
ALLOBJ := $(ALLOBJ:.f90=.o)
ALLOBJ := $(ALLOBJ:.F90=.o)
ALLOBJ := $(ALLOBJ:.f=.o)

EXE = $(ALLOBJ:.o=)

all: $(EXE)

$(EXE) : $(LIBCOREBLAS) $(LIBPLASMA) $(QUARKDIR)/libquark.a 

clean:
	rm -f $(ALLOBJ) 

cleanall: clean
	rm -f $(EXE)

%.o : %.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

%.o : %.F90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.f.o:
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

example_% : example_%.o
	$(LOADER) $(LDFLAGS) $< -o $@ $(LIB)

.PHONY: clean cleanall
