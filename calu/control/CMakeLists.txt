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

# THE LIBRARY IS BUILT IN THE COMPUTE SUBDIRECTORY

ADD_CUSTOM_TARGET( control_include ALL DEPENDS
  compute_d.h compute_c.h compute_s.h )

