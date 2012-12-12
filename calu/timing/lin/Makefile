###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.5.0
# @author Piotr Luszczek
# @date 2010-11-15
#
###

PLASMA_DIR = ../..
include ../../Makefile.internal

.SUFFIXES:
.SUFFIXES: .f .f90 .o

ALINTIM = atimin.o atimck.o icopy.o xlaenv.o

DZINTIM = dlaord.o dprtbl.o dprtb2.o dprtb3.o dprtb4.o dprtb5.o \
          dprtls.o dmflop.o dopaux.o dopbl2.o dopbl3.o dopla.o

DLINTIM = gnufor2.o drawer.o dtimaa.o dtimmg.o \
          dtimmv.o dtimmm.o dtimb2.o dtimb3.o dtimge.o \
          dtimls.o dtimpo.o \
          dtimqr.o dtimlq.o \
          dqrt13.o dqrt15.o dlinpk.o deispk.o

all:  double 

double: xlintimd

INC        = -I../../include
LIBDIR     = -L../../lib

xlintimd : $(ALINTIM) $(DZINTIM) $(DLINTIM)
	$(FC) $(LDFLAGS) $(ALINTIM) $(DLINTIM) $(DZINTIM) $(LIBDIR) $(LIB) -o $@

$(ALINTIM): $(FRC)
$(DZINTIM): $(FRC)
$(DLINTIM): $(FRC)

FRC:
	@FRC=$(FRC)

clean:
	rm -f *.o *.mod *.txt *.ps

.f90.o:
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

dtimaa.o: dtimaa.f
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

.f.o:
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

.o:
	$(FC) $(LDFLAGS) $(INC) $< -o $@ $(LIBDIR) $(LIB) 
