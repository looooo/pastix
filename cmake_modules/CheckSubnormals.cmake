# cmake modules setup
cmake_minimum_required (VERSION 3.5)
include(CheckCSourceCompiles)
#
# Subnormals
#
check_c_source_compiles("
#include <xmmintrin.h>
int main(void){
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    return 0;
}"
HAVE_FTZ_MACROS
)

check_c_source_compiles("
#include <pmmintrin.h>
int main(void){
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    return 0;
}"
HAVE_DAZ_MACROS
)

if (NOT HAVE_FTZ_MACROS AND NOT HAVE_DAZ_MACROS)
    check_c_source_compiles("
    #include <immintrin.h>
    int main(void){
    _mm_setcsr(_mm_getcsr() & ~0x8040);
    return 0;
    }"
    HAVE_MM_SETCSR
    )
endif (NOT HAVE_FTZ_MACROS AND NOT HAVE_DAZ_MACROS)