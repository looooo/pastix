###
#
# @file Makefile
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.6
# @author Julie Langou
# @date 2010-11-15
#
###

.PHONY: routine htmlbrowsing

all: routine htmlbrowsing

reference_guide:
	(cd ../../ && make generate && cd docs/src)
	./get_comments.sh > reference_guide.txt
	asciidoc -a toc reference_guide.txt
	mv -f reference_guide.html ../html
	sed -i -e "s/echo \"\*\*\*\*\"//g" reference_guide.txt
#	Transform optional title in title for the pdf
	sed -i -e "s/^\./=== /g" reference_guide.txt
	asciidoc -a docinfo -b docbook reference_guide.txt
	dblatex  reference_guide.xml
# 	get a letter format of the A4 pdf
#	pdfcrop --margins 15 reference_guide.pdf ../pdf/routineLetter.pdf
	mv -f reference_guide.pdf ../pdf
	rm -f reference_guide.txt
	rm -f reference_guide.xml

htmlbrowsing:
	(cd ../../ && make generate && cd docs/src)
	rm -f ../html/htmlbrowsing/*
	cp ../../src/*.c ../../cblas/*.c ../html/htmlbrowsing
	(cd ../html/htmlbrowsing ; source-highlight  -n --gen-references=inline   --ctags='ctags --excmd=n --tag-relative=yes --file-scope=no' *.c ; rm -f *.c tags ; cd -)
	cp ../../core_blas/*.f ../../core_lapack/*.f ../html/htmlbrowsing
	(cd ../html/htmlbrowsing ; source-highlight --src-lang=fortran -n --gen-references=inline   --ctags='ctags --excmd=n --tag-relative=yes --file-scope=no' *.f ; rm -f *.f tags ; cd -)
	./add_FortranBrowsing.sh
