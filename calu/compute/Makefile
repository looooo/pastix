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

# Overwritten in make.inc
PLASMA_DIR = ..
include ../Makefile.internal

INC := -I$(PLASMA_DIR)/control/ $(INC)

ZHDR =  zgecfi2.h

ZSRC =  \
	zgels.c  zgelqf.c zgelqs.c zgeqrf.c zgeqrs.c		\
	zgesv_incpiv.c  zgetrf_incpiv.c zgetrs_incpiv.c		\
	zgesv.c  zgetrf.c zgetrs.c zlaswp.c zlaswpc.c		\
	zposv.c  zpotrf.c zpotrs.c				\
	ztrsm.c  zgemm.c  zsymm.c zsyrk.c ztrsmpl.c ztrsmrv.c	\
	zpotri.c ztrtri.c zlauum.c ztrmm.c zsyr2k.c		\
	zgetri.c						\
	zlange.c zlansy.c zlaset.c zlacpy.c			\
	zunmqr.c zunmlq.c zungqr.c zunglq.c			\
	zhemm.c  zherk.c  zlanhe.c zher2k.c			\
	zcposv.c zcgesv.c zcungesv.c				\
	zgecfi2.c zgecfi.c zgetmi.c				\
	zplrnt.c zplghe.c zplgsy.c				\
	zhegv.c zhegvd.c  zheev.c zheevd.c zhegst.c zhetrd.c	\
	zgesvd.c zgebrd.c					\
	zgecon.c zpocon.c					\
	zgeqp3.c						\
	\
	pzgelqf.c pzgeqrf.c pzpotrf.c pztrtri.c pzlauum.c		\
	pzgelqfrh.c pzgeqrfrh.c						\
	pztrsm.c  pzgemm.c  pzsymm.c  pzsyrk.c  pztrsmpl.c pztrsmrv.c	\
	pzlange.c pzlansy.c pzsyr2k.c					\
	pzlacpy.c pzgeadd.c  pzgetrf_incpiv.c pztrmm.c			\
	pzunmqr.c pzunmlq.c pzungqr.c pzunglq.c				\
	pzunmqrrh.c pzungqrrh.c						\
	pzunmlqrh.c pzunglqrh.c						\
	pzlanhe.c pzhemm.c pzherk.c pzher2k.c				\
	pzlag2c.c pzlaset.c pzlaset2.c					\
	pzpack.c  pzshift.c pzgetmi2.c					\
	pzplrnt.c pzplghe.c pzplgsy.c					\
	pzgetrf_reclap.c pzgetrf_rectil.c pzlaswp.c pzlaswpc.c		\
	pzbarrier.c							\
	pzgebrd_ge2tb.c pzgebrd_tb2bd.c					\
	pzhetrd_he2hb.c pzhetrd_hb2st.c					\
	pzhegst.c							\
	pzlarft_blgtrd.c pzunmqr_blgtrd.c				\
	pzhbcpy_t2bl.c							\
	pzgeqp3.c							\

ZSRC =  \
        zgesv_incpiv.c  zgetrf_incpiv.c zgetrs_incpiv.c		\
        zgesv.c  zgetrf.c zgetrs.c zlaswp.c zlaswpc.c		\
	ztrsm.c  zgemm.c  zsymm.c zsyrk.c ztrsmpl.c ztrsmrv.c  	\
	zlange.c zlansy.c zlaset.c zlacpy.c			\
	zhemm.c  zherk.c  zlanhe.c zher2k.c   			\
        zgecfi2.c zgecfi.c zgetmi.c                     	\
	zplrnt.c zplghe.c zplgsy.c				\
	\
        pztrsm.c  pzgemm.c  pzsymm.c  pzsyrk.c  pztrsmpl.c pztrsmrv.c 	\
	pzlange.c pzlansy.c pzsyr2k.c					\
	pzlacpy.c pzgeadd.c  pzgetrf_incpiv.c pztrmm.c			\
	pzlanhe.c pzhemm.c pzherk.c pzher2k.c				\
	pzlag2c.c pzlaset.c pzlaset2.c                          	\
	pzpack.c  pzshift.c pzgetmi2.c					\
	pzplrnt.c pzplghe.c pzplgsy.c                           	\
	pzgetrf_reclap.c pzgetrf_rectil.c pzlaswp.c pzlaswpc.c		\
	pzbarrier.c                                             	\
	pzgetrf_nopiv.c zgetrf_nopiv.c                                  \
        pzgetrf_tntpiv.c

-include Makefile.src

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.c=.o)

all: $(LIBPLASMA)

$(LIBPLASMA): $(HEADERS) $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean:
	rm -f *.o *~

cleanall: clean
	rm -f *.a

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

zgecfi2.c zgecfi.c : $(ZHDR)
cgecfi2.c cgecfi.c : $(CHDR)
dgecfi2.c dgecfi.c : $(DHDR)
sgecfi2.c sgecfi.c : $(SHDR)

.PHONY: clean cleanall
