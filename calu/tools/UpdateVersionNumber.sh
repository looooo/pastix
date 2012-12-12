#!/bin/sh

major=2
minor=5
micro=0

find -name "*.[ch]" | xargs sed -i 's/@version 2.4.6/@version 2.5.0/'
find -name "CMakeLists.txt" | xargs sed -i 's/@version 2.4.6/@version 2.5.0/'
find -name "Makefile" | xargs sed -i 's/@version 2.4.6/@version 2.5.0/'
sed -i 's/@version 2.4.6/@version 2.5.0/' Makefile.gen Makefile.internal Makefile.tau
sed -i 's/@version 2.4.6/@version 2.5.0/' makes/make.inc.*
sed -i 's/@version 2.4.6/@version 2.5.0/' make.inc.example
sed -i 's/@version 2.4.6/@version 2.5.0/' docs/latex/contributors_guide/comments.tex
sed -i 's/@version 2.4.6/@version 2.5.0/' makes/cmake32.bat makes/cmake64.bat
sed -i 's/@version 2.4.6/@version 2.5.0/' tools/convert2eztrace.pl
sed -i 's/@version 2.4.6/@version 2.5.0/' tools/genf77interface.pl
sed -i 's/@version 2.4.6/@version 2.5.0/' tools/genf90interface.pl
sed -i 's/version 2.4.6/version 2.5.0/' examples/*.f
sed -i 's/Version: 2.4.6/Version: 2.5.0/' lib/pkgconfig/plasma.pc.in

sed -i 's/VERSION \?= 2.4.6/VERSION ?= 2.5.0/' docs/doxygen/Makefile

sed -i 's/PLASMA_VERSION_MAJOR[ ]*[0-9]/PLASMA_VERSION_MAJOR 2/' include/plasma.h
sed -i 's/PLASMA_VERSION_MINOR[ ]*[0-9]/PLASMA_VERSION_MINOR 4/' include/plasma.h
sed -i 's/PLASMA_VERSION_MICRO[ ]*[0-9]/PLASMA_VERSION_MICRO 6/' include/plasma.h

sed -i 's/PLASMA_VERSION_MAJOR[ ]*"[0-9]"/PLASMA_VERSION_MAJOR "2"/'  CMakeLists.txt
sed -i 's/PLASMA_VERSION_MINOR[ ]*"[0-9]"/PLASMA_VERSION_MINOR "5"/'  CMakeLists.txt
sed -i 's/PLASMA_VERSION_PATCH[ ]*"[0-9]"/PLASMA_VERSION_PATCH "0"/'  CMakeLists.txt


sed -i 's/version 2.4.6/version 2.5.0/' ../plasma-installer/script/*.py
sed -i 's/version 2.4.6/version 2.5.0/' ../plasma-installer/setup.py
sed -i 's/"2.4.6"/"2.4.6", "2.5.0"/' ../plasma-installer/script/framework.py
