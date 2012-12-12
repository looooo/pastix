###
#
# @file CMakeLists.txt
#
# PLASMA is a software package provided by Univ. of Tennessee,
# Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.4.6
# @author Asim YarKhan
# @date 2010-11-15
#
###

INCLUDE( ${CMAKE_SOURCE_DIR}/GenerateZDCS.cmake OPTIONAL )

set( NAMESC
        gelqf
	gels
	geqrs
	getrs
	potrf
	trsm
        gelqs
	geqrf
	gesv
	posv
	potrs
	trsmpl
##	unmqr
##	unmlq
)

set( NAMESF
	posv_f
	gesv_f
	gels_f
)

foreach( name ${NAMESC} )
  set( zCEXE ${zCEXE} example_z${name} )
  set( cCEXE ${cCEXE} example_c${name} )
  set( dCEXE ${dCEXE} example_d${name} )
  set( sCEXE ${sCEXE} example_s${name} )
endforeach()
set( zCEXE ${zCEXE} example_zunmqr example_zunmlq )
set( cCEXE ${cCEXE} example_cunmqr example_cunmlq )
set( dCEXE ${dCEXE} example_dormqr example_dormlq )
set( sCEXE ${sCEXE} example_sormqr example_sormlq )
set( CEXE ${zCEXE} ${cCEXE} ${dCEXE} ${sCEXE} )

foreach( name ${NAMESF} )
  set( zFEXE ${zFEXE} example_z${name} )
  set( cFEXE ${cFEXE} example_c${name} )
  set( dFEXE ${dFEXE} example_d${name} )
  set( sFEXE ${sFEXE} example_s${name} )
endforeach()
set( FEXE ${zFEXE} ${cFEXE} ${dFEXE} ${sFEXE} )

FOREACH( name ${CEXE} )
  set( CEXE_SRC ${CEXE_SRC} ${name}.c )
  ADD_EXECUTABLE( ${name} ${name}.c )
  PLASMA_SETUP_C_MAIN( ${name} )
ENDFOREACH()

FOREACH( name ${FEXE} )
  set( FEXE_SRC ${FEXE_SRC} ${name}.f )
  ADD_EXECUTABLE( ${name} ${name}.f )
  PLASMA_SETUP_FORTRAN_MAIN( ${name} )
ENDFOREACH()

ADD_CUSTOM_TARGET( examples_SRC ALL DEPENDS ${CEXE_SRC} ${FEXE_SRC} )

IF(WIN32 AND NOT UNIX)
  INSTALL( FILES Makefile.nmake CMakeLists.txt DESTINATION examples )
ELSE()
  INSTALL( FILES CMakeLists.txt DESTINATION examples )
ENDIF()

FOREACH( name ${CEXE_SRC} ${FEXE_SRC} )
  INSTALL( FILES ${name} DESTINATION examples )
ENDFOREACH()

#INSTALL( TARGETS ${CEXE} ${FEXE} DESTINATION examples )

