###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.6
# @author Bilel Hadri
# @author Hatem Ltaief
# @date 2010-11-15
#
###

PLASMA_DIR = ..
include ../Makefile.internal

ZHDR=	testing_zmain.h

ZSRC=	testing_zmain.c 	\
	testing_zgesv.c 	\
	testing_zgetri.c 	\
	testing_zgesvd.c 	\
	testing_zgels.c 	\
	testing_zposv.c 	\
	testing_zpotri.c 	\
	testing_zgemm.c 	\
	testing_zheev.c 	\
	testing_zheevd.c 	\
	testing_zhegst.c 	\
	testing_zhegv.c 	\
	testing_zhegvd.c 	\
	testing_zsyrk.c 	\
	testing_zsyr2k.c 	\
	testing_zsymm.c 	\
	testing_ztrsm.c 	\
	testing_ztrmm.c 	\
	testing_zlange.c	\
	testing_zherk.c		\
	testing_zher2k.c	\
	testing_zhemm.c		\
	testing_zgecfi.c	\
	testing_zgetmi.c	\
	testing_zcposv.c 	\
	testing_zcgesv.c 	\
	testing_zcungesv.c  \
	testing_zpemv.c 

-include Makefile.src

ALLHDR= $(ZHDR) $(CHDR) $(DHDR) $(SHDR)

ZOBJ = $(ZSRC:.c=.o)
COBJ = $(CSRC:.c=.o)
DOBJ = $(DSRC:.c=.o)
SOBJ = $(SSRC:.c=.o)

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.c=.o)

EXE = ztesting ctesting dtesting stesting

all: $(EXE)

$(EXE) : $(LIBCOREBLAS) $(LIBPLASMA) $(QUARKDIR)/libquark.a 

clean:
	rm -f $(ALLOBJ) *~

cleanall: clean
	rm -f $(EXE)

testing_z%.o: testing_z%.c testing_zmain.h
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

testing_c%.o: testing_c%.c testing_cmain.h
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

testing_d%.o: testing_d%.c testing_dmain.h
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

testing_s%.o: testing_s%.c testing_smain.h
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.f.o:
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

ztesting : testing_zmain.h $(ZOBJ)
	$(LOADER) $(LDFLAGS) $(ZOBJ) -o $@ $(LIB)

ctesting : testing_cmain.h $(COBJ)
	$(LOADER) $(LDFLAGS) $(COBJ) -o $@ $(LIB)

dtesting : testing_dmain.h $(DOBJ)
	$(LOADER) $(LDFLAGS) $(DOBJ) -o $@ $(LIB)

stesting : testing_smain.h $(SOBJ)
	$(LOADER) $(LDFLAGS) $(SOBJ) -o $@ $(LIB)

.PHONY: clean cleanall
