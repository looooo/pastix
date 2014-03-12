###
#
# -- Inria
# -- (C) Copyright 2012
#
# This software is a computer program whose purpose is to process
# Matrices Over Runtime Systems @ Exascale (MORSE). More information
# can be found on the following website: http://www.inria.fr/en/teams/morse.
# 
# This software is governed by the CeCILL-C license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL-C
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
# 
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
# 
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
# 
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL-C license and that you accept its terms.
#
###
#
#  @file ParseArguments.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver. 
# 
#  @version 2.0.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @date 13-07-2012
#   
###

MACRO(PARSE_ARGUMENTS prefix arg_names option_names)
    set(DEFAULT_ARGS)
    foreach(arg_name ${arg_names})    
        set(${prefix}_${arg_name})
    endforeach(arg_name)
    foreach(option ${option_names})
        set(${prefix}_${option} FALSE)
    endforeach(option)

    set(current_arg_name DEFAULT_ARGS)
    set(current_arg_list)
    foreach(arg ${ARGN})        
        set(larg_names ${arg_names})    
        list(FIND larg_names "${arg}" is_arg_name)           
        if(is_arg_name GREATER -1)
            set(${prefix}_${current_arg_name} ${current_arg_list})
            set(current_arg_name ${arg})
            set(current_arg_list)
        else(is_arg_name GREATER -1)
            set(loption_names ${option_names})    
            list(FIND loption_names "${arg}" is_option)        
            if(is_option GREATER -1)
                set(${prefix}_${arg} TRUE)
            else(is_option GREATER -1)
                set(current_arg_list ${current_arg_list} ${arg})
            endif(is_option GREATER -1)
        endif(is_arg_name GREATER -1)
    endforeach(arg)
    set(${prefix}_${current_arg_name} ${current_arg_list})

ENDMACRO(PARSE_ARGUMENTS)

MACRO(CAR var)
    set(${var} ${ARGV1})
ENDMACRO(CAR)

MACRO(CDR var junk)
    set(${var} ${ARGN})
ENDMACRO(CDR)

##
## @end file ParseArguments.cmake
##
