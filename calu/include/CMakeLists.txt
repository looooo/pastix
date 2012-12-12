###
#
# @file CMakeLists.txt
#
# PLASMA is a software package provided by Univ. of Tennessee,
# Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.5.0
# @author Asim YarKhan
# @date 2010-11-15
#
###

INCLUDE( ${CMAKE_SOURCE_DIR}/GenerateZDCS.cmake OPTIONAL )

ADD_CUSTOM_TARGET( include ALL DEPENDS 
   plasma_s.h plasma_d.h plasma_c.h 
   plasma_ds.h 
   core_sblas.h core_dblas.h core_cblas.h 
   core_dsblas.h )

INSTALL(
  FILES cblas.h core_cblas.h core_dsblas.h core_zblas.h plasma_c.h
    plasma_ds.h plasma.h plasma_z.h core_blas.h
    core_dblas.h core_sblas.h core_zcblas.h plasma_d.h plasmaf.h
    plasma_s.h plasma_zc.h CMakeLists.txt
  DESTINATION include )
