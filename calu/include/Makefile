PLASMA_DIR = ..
include ../Makefile.internal

ZHDR =  plasma_z.h plasma_zc.h	\
	core_zblas.h core_zcblas.h
CHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p c)
DHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p "d ds")
SHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p s)

ALLHDR= cblas.h plasmaf.h	\
	core_blas.h plasma.h	\
	$(ZHDR) $(CHDR) $(DHDR) $(SHDR)


all: $(ALLHDR)

clean:
	rm -f *.o *~

cleanall: clean
	rm -f *.a
