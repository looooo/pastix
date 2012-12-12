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

.PHONY: users_guide contributors_guide plasma_tau

all: users_guide contributors_guide plasma_tau

users_guide:
	( cd users_guide && pdflatex users_guide.tex && bibtex users_guide && pdflatex users_guide.tex && pdflatex users_guide.tex )
	mv users_guide/users_guide.pdf ../pdf/users_guide.pdf
	( cd users_guide && rm -f *.aux *.log *.out *.toc *.bbl *.blg )

contributors_guide:
	( cd contributors_guide && pdflatex contributors_guide.tex && pdflatex contributors_guide.tex )
	mv contributors_guide/contributors_guide.pdf ../pdf/contributors_guide.pdf
	( cd contributors_guide && rm -f *.aux *.log *.out *.toc )

plasma_tau:
	( cd plasma_tau && pdflatex plasma_tau.tex && pdflatex plasma_tau.tex )
	mv plasma_tau/plasma_tau.pdf ../pdf/plasma_tau.pdf
	( cd plasma_tau && rm -f *.aux *.log *.out *.toc )

clean:
	rm -f ../pdf/users_guide.pdf
	rm -f ../pdf/contributors_guide.pdf
	rm -f ../pdf/plasma_tau.pdf
