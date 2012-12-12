###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.5.0
# @author Jakub Kurzak
# @date 2010-11-15
#
###

PLASMA_DIR = ..
include ../Makefile.internal

ZHDR=	compute_z.h

ZSRC=	plasma_zf77.c workspace_z.c 	\
	pztile.c ztile.c

ZF90=	plasma_zf90.F90 		\
	plasma_zcf90.F90 		\
	plasma_zf90_wrappers.F90

-include Makefile.src

ALLHDR= allocate.h   async.h   auxiliary.h	\
	common.h     context.h control.h	\
        descriptor.h global.h  tile.h		\
        workspace.h  primes.h gkkleader.h       \
	$(ZHDR) $(CHDR) $(DHDR) $(SHDR)

ALLSRC= allocate.c async.c    auxiliary.c	\
	context.c  control.c  descriptor.c 	\
	tile.c     plasmaos.c workspace.c	\
        primes.c   gkkleader.c                  \
	plasma_f77.c plasma_mf77.c		\
	$(ZSRC) $(CSRC) $(DSRC) $(SSRC)

ifeq (${PLASMA_F90}, 1)
ALLSRC:= $(ALLSRC)	\
	 plasma_f90.f90	\
	 $(ZF90) $(CF90) $(DF90) $(SF90)

F90MOD= plasma_zc.mod 	\
	plasma_ds.mod 	\
	plasma_z.mod 	\
	plasma_c.mod	\
	plasma_d.mod 	\
	plasma_s.mod	\
	plasma.mod

ALLMOD= ${PLASMA_DIR}/include/plasma_zc.mod 	\
	${PLASMA_DIR}/include/plasma_ds.mod	\
	${PLASMA_DIR}/include/plasma_z.mod 	\
	${PLASMA_DIR}/include/plasma_c.mod	\
	${PLASMA_DIR}/include/plasma_d.mod 	\
	${PLASMA_DIR}/include/plasma_s.mod	\
	${PLASMA_DIR}/include/plasma.mod

endif

ALLOBJ  = $(ALLSRC:.c=.o)
ALLOBJ := $(ALLOBJ:.f90=.o)
ALLOBJ := $(ALLOBJ:.F90=.o)

all: $(LIBPLASMA)

ifeq (${PLASMA_F90}, 1)
$(LIBPLASMA) : $(ALLMOD)
endif

$(LIBPLASMA): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean:
	rm -f *.o *~ *.mod

cleanall: clean
	rm -f *.a $(ALLMOD)

plasmaos.o : plasmaos.c plasmaos-hwloc.c 

plasma_zc.mod : plasma_zcf90.o
plasma_ds.mod : plasma_dsf90.o
plasma_z.mod : plasma_zf90.o
plasma_c.mod : plasma_cf90.o
plasma_d.mod : plasma_df90.o
plasma_s.mod : plasma_sf90.o
plasma.mod   : plasma_f90.o 
plasma_f90.o : plasma_zf90.o plasma_cf90.o plasma_df90.o plasma_sf90.o plasma_dsf90.o plasma_zcf90.o
plasma_zf90_wrappers.o plasma_cf90_wrappers.o plasma_df90_wrappers.o plasma_sf90_wrappers.o : ${PLASMA_DIR}/include/plasma.mod

${PLASMA_DIR}/include/%.mod : %.mod
	cp $< $@

.c.o: $(ALLHDR)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

%.o : %.F90
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

.f.o:
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

# Headers dependencies
$(ALLOBJ) : $(ALLHDR)

.PHONY: clean cleanall
