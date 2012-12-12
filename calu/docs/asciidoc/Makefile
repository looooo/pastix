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

.PHONY: readme release_notes todo cmake_build_notes installguide

all: readme release_notes todo cmake_build_notes installation_guide

readme:
	asciidoc -n -a toc README.txt
	mv -f README.html ../html/README.html
	asciidoc -n -a toc -b html4 -o - README.txt | lynx -justify -dump -stdin > ../../README

release_notes:
	asciidoc -a toc ReleaseNotes.txt
	mv -f ReleaseNotes.html ../html/ReleaseNotes.html
	asciidoc -a toc -b html4 -o - ReleaseNotes.txt | lynx -justify -dump -stdin > ../../ReleaseNotes

todo:
	asciidoc -a toc ToDo.txt
	mv -f ToDo.html ../html/ToDo.html
	asciidoc -a toc -b html4 -o - ToDo.txt | lynx -justify -dump -stdin > ../../ToDo

cmake_build_notes:
	asciidoc -a toc CMakeBuildNotes.txt
	mv -f CMakeBuildNotes.html ../html/CMakeBuildNotes.html
	asciidoc -a toc -b html4 -o - CMakeBuildNotes.txt | lynx -width=110 -justify -dump -stdin > ../../CMakeBuildNotes

installation_guide:
	asciidoc -a toc InstallationGuide.txt
	mv -f InstallationGuide.html ../html/InstallationGuide.html
	asciidoc -a toc -b html4 -o - InstallationGuide.txt | lynx -justify -dump -stdin > ../../InstallationGuide

clean:
	rm -f ../html/README.html
	rm -f ../html/ReleaseNotes.html
	rm -f ../html/ToDo.html
	rm -f ../html/CMakeBuildNotes.html
	rm -f ../html/InstallationGuide.html
	rm -f ../../README
	rm -f ../../ReleaseNotes
	rm -f ../../ToDo
	rm -f ../../CMakeBuildNotes
	rm -f ../../InstallationGuide
