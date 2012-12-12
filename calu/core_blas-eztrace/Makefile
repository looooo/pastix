###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.6
# @date 2010-11-15
#
###

PLASMA_DIR = ..
include ../Makefile.internal

# Add the headers from src directory
INC :=  -I../control/ $(INC) \
	-I$(EZT_DIR)/include -I$(GTG_DIR)/include -I$(FXT_DIR)/include

ZSRC=	coreblas_z.c

-include Makefile.src

ALLSRC = coreblas.c \
	 $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.c=.o)

all: $(LIBEZT_COREBLAS) $(LIBEZT_CONVERT)

$(LIBEZT_COREBLAS) : $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

$(LIBEZT_CONVERT): eztrace_convert_coreblas.o
	$(CC) -shared $(CFLAGS) $(INC) $^ -o $@

clean:
	rm -f *.o *~

cleanall: clean
	rm -f $(LIBEZT_COREBLAS) $(LIBEZT_CONVERT)

eztrace_convert_coreblas.o :  eztrace_convert_coreblas.c
	$(CC) -fPIC $(CFLAGS) $(INC) -c $< -o $@

.c.o :
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.PHONY: all clean cleanall
