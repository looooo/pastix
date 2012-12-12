###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.5.0
# @author Mathieu Faverge
# @date 2010-11-15
#
###

PLASMA_DIR = ..
include ../Makefile.internal

CFLAGS  += -I../control/ $(INC)
LDFLAGS += $(LIB)

ZHDR=	zauxiliary.h zauxiliary.c

ZSRC=	time_zgemm.c 		\
	time_ztrsm.c 		\
	time_zcposv.c 		\
	time_zposv.c 		\
	time_zpotrf.c 		\
	time_zcgesv.c 		\
	time_zgesv.c 		\
	time_zgetrf.c 		\
	time_zgesv_incpiv.c 	\
	time_zgetrf_incpiv.c	\
	time_zgels.c 		\
	time_zgeqrf.c 		\
	time_zgemm_tile.c 	\
	time_zcposv_tile.c 	\
	time_zposv_tile.c 	\
	time_zpotrf_tile.c 	\
	time_zpotri_tile.c	\
	time_zcgesv_tile.c 	\
	time_zgesv_tile.c 	\
	time_zgetrf_tile.c 	\
	time_zgetri_tile.c 	\
	time_zgesv_incpiv_tile.c 	\
	time_zgetrf_incpiv_tile.c 	\
	time_zgels_tile.c 	\
	time_zgeqrf_tile.c      \
        time_zgecfi.c           \
        time_zlapack2tile.c	\
	\
	time_zheev_tile.c	\
	time_zheevd_tile.c	\
	time_zhegv_tile.c	\
	time_zgesvd_tile.c	\
	time_zgetrf_reclap.c	\
	time_zgetrf_rectil.c	\
	time_zgeqp3.c		\

ZSRC=	time_zgesv.c 		\
	time_zgetrf.c 		\
	time_zgesv_incpiv.c 	\
	time_zgetrf_incpiv.c	\
	time_zgesv_tile.c 	\
	time_zgetrf_tile.c 	\
	time_zgesv_incpiv_tile.c 	\
	time_zgetrf_incpiv_tile.c 	\
	time_zgetrf_reclap.c	\
	time_zgetrf_rectil.c

-include Makefile.src

#timing.c is in headers because it's a template for every timing routine
ALLHDR= auxiliary.h auxiliary.c timing.h timing.c \
	$(ZHDR) $(CHDR) $(DHDR) $(SHDR)

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.c=.o)
EXE = $(ALLOBJ:.o=)

all: $(EXE)

$(EXE) : $(LIBCOREBLAS) $(LIBPLASMA) $(QUARKDIR)/libquark.a 

clean:
	rm -f *.o *~

cleanall: clean
	rm -f $(EXE)

%.c: $(ALLHDR)

time_%.o : time_%.c $(ALLHDR)
	$(CC) $(CFLAGS) -c $< -o $@

%auxiliary.o : %auxiliary.c $(ALLHDR)
	$(CC) $(CFLAGS) -c $< -o $@

time_s% : time_s%.o sauxiliary.o
	$(LOADER) $^ -o $@ $(LDFLAGS)
time_d% : time_d%.o dauxiliary.o
	$(LOADER) $^ -o $@ $(LDFLAGS)
time_c% : time_c%.o cauxiliary.o
	$(LOADER) $^ -o $@ $(LDFLAGS)
time_z% : time_z%.o zauxiliary.o
	$(LOADER) $^ -o $@ $(LDFLAGS)

time_zlapack2tile time_clapack2tile time_dlapack2tile time_slapack2tile : auxiliary.o
time_zgecfi time_cgecfi time_dgecfi time_sgecfi : auxiliary.o

png: time_dgesv time_zgesv
	./time_dgesv --gnuplot-png | tee timeplot.gpt ; gnuplot timeplot.gpt

.PHONY: all clean

.SECONDARY: $(ALLOBJ)
