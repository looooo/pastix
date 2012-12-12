REM ###
REM #
REM # @file cmake32.bat
REM #
REM #  PLASMA is a software package provided by Univ. of Tennessee,
REM #  Univ. of California Berkeley and Univ. of Colorado Denver
REM #
REM # @version 2.4.6
REM # @author Piotr Luszczek
REM # @date 2010-11-15
REM #
REM ###

REM Build using CMake on a Windows 32 machine

REM Put CMake on the path
SET PATH=c:/Program Files (x86)/CMake 2.6/bin;%PATH%
REM rename GenerateZDCS.cmake.optional
IF EXIST GenerateZDCS.cmake DEL GenerateZDCS.cmake
IF EXIST GenerateZDCS.cmake.optional RENAME GenerateZDCS.cmake.optional GenerateZDCS.cmake
REM Setup and build 32 bit Intel Fortran env
CALL "c:/Program Files (x86)/Intel/Compiler/Fortran/10.1.025/IA32/Bin/ifortvars.bat"
REM Setup and build 32 bit Intel C++ env
CALL "c:/Program Files (x86)/Intel/Compiler/C++/10.1.025/IA32/Bin/iclvars.bat"
REM Setup and build 32 bit Intel MKL env
CALL "c:/Program Files/Intel/MKL/10.0.4.023/tools/environment/mklvars32.bat"
REM Cleanup Cmake files
REM del CMakeCache.txt
REM Make build_dir if not present and cd to build_dir
IF EXIST build_dir (
DEL /F /Q build_dir
) ELSE (
MKDIR build_dir
)
CD build_dir
REM Call CMake - remember to fix the path to the source
cmake -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_C_COMPILER=icl -DCMAKE_CXX_COMPILER=icl -G "NMake Makefiles" D:/fike/buildbot/plasma/mordor_win32_intel_mkl/build
nmake package

REM Do the usual testing ...
