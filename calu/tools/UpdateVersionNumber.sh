#!/bin/sh

major=2
minor=4
micro=6

find -name "*.[ch]" | xargs sed -i 's/@version 2.4.5/@version 2.4.6/'
find -name "CMakeLists.txt" | xargs sed -i 's/@version 2.4.5/@version 2.4.6/'
find -name "Makefile" | xargs sed -i 's/@version 2.4.5/@version 2.4.6/'
sed -i 's/@version 2.4.5/@version 2.4.6/' Makefile.gen Makefile.internal Makefile.tau 
sed -i 's/@version 2.4.5/@version 2.4.6/' makes/make.inc.*
sed -i 's/@version 2.4.5/@version 2.4.6/' make.inc.example                            
sed -i 's/@version 2.4.5/@version 2.4.6/' docs/latex/contributors_guide/comments.tex
sed -i 's/@version 2.4.5/@version 2.4.6/' makes/cmake32.bat makes/cmake64.bat       
sed -i 's/@version 2.4.5/@version 2.4.6/' tools/convert2eztrace.pl           
sed -i 's/@version 2.4.5/@version 2.4.6/' tools/genf77interface.pl
sed -i 's/@version 2.4.5/@version 2.4.6/' tools/genf90interface.pl
sed -i 's/version 2.4.5/version 2.4.6/' examples/*.f
sed -i 's/Version: 2.4.5/Version: 2.4.6/' lib/pkgconfig/plasma.pc.in

sed -i 's/VERSION \?= 2.4.5/VERSION ?= 2.4.6/' docs/doxygen/Makefile

sed -i 's/PLASMA_VERSION_MAJOR[ ]*[0-9]/PLASMA_VERSION_MAJOR 2/' include/plasma.h
sed -i 's/PLASMA_VERSION_MINOR[ ]*[0-9]/PLASMA_VERSION_MINOR 4/' include/plasma.h
sed -i 's/PLASMA_VERSION_MICRO[ ]*[0-9]/PLASMA_VERSION_MICRO 6/' include/plasma.h

sed -i 's/PLASMA_VERSION_MAJOR[ ]*"[0-9]"/PLASMA_VERSION_MAJOR "2"/'  CMakeLists.txt
sed -i 's/PLASMA_VERSION_MINOR[ ]*"[0-9]"/PLASMA_VERSION_MINOR "4"/'  CMakeLists.txt
sed -i 's/PLASMA_VERSION_PATCH[ ]*"[0-9]"/PLASMA_VERSION_PATCH "6"/'  CMakeLists.txt


sed -i 's/version 2.4.5/version 2.4.6/' ../plasma-installer/script/*.py
sed -i 's/version 2.4.5/version 2.4.6/' ../plasma-installer/setup.py
sed -i 's/"2.4.5"/"2.4.6", "2.4.5"/' ../plasma-installer/script/framework.py
