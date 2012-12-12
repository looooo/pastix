###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.5.0
# @author Wesley Alvaro
# @date 2010-11-15
#
###

DOXYGEN ?= doxygen
VERSION ?= 2.4.6
PLASMA_DIR = ../..
DOCS_DIR = $(PLASMA_DIR)/docs
.PHONY: clean realclean doxygen doxygen-usr doxygen-dev install-latex install-html install all again

all: doxygen
again: clean doxygen

install: install-html install-latex

install-html: doxygen
	( test ! -d ../html/doxygen && mkdir ../html/doxygen ) || true
	cp -rf out/html/* $(DOCS_DIR)/html/doxygen/.

install-latex: doxygen
	( test ! -d $(DOCS_DIR)/latex/doxygen && mkdir $(DOCS_DIR)/latex/doxygen ) || true
	( cd out/latex && $(MAKE) )
	cp out/latex/*.pdf $(DOCS_DIR)/latex/doxygen/.

doxygen: out

out: plasma.dox groups-usr.dox groups-dev.dox
	$(DOXYGEN) plasma.dox

plasma.dox:
	if	test -f $(PLASMA_DIR)/include/Makefile; \
		then make plasma-dev.dox; \
		else make plasma-usr.dox; \
	fi

plasma-usr.dox: conf-usr.dox
	sed s/!VERSION!/$(VERSION)/g conf-usr.dox > plasma-usr.dox
	cp plasma-usr.dox plasma.dox

plasma-dev.dox: conf-dev.dox
	sed s/!VERSION!/$(VERSION)/g conf-dev.dox > plasma-dev.dox
	cp plasma-dev.dox plasma.dox

realclean: clean
	rm -rf $(DOCS_DIR)/html/doxygen
	rm -rf $(DOCS_DIR)/latex/doxygen

clean:
	rm -rf plasma-usr.dox
	rm -rf plasma-dev.dox
	rm -rf plasma.dox
	rm -rf out

