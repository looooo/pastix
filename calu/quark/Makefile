###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.5.0
# @author Asim YarKhan
# @date 2010-11-15
#
###

PLASMA_DIR = ..
include ../Makefile.internal

all: 	libquark.a

lib: 	libquark.a

clean:
	rm -f *.o *~

cleanall: clean
	rm -f *.a

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

libquark.a: icl_hash.o icl_list.o quarkos.o quark.o
	$(ARCH) $(ARCHFLAGS) $@ $^
	$(RANLIB) $@

.PHONY: clean

