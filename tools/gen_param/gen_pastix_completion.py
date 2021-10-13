"""
 @file gen_pastix_completion.py

 Generate the pastix_completion.bash script that allows auto-completion
 of pastix execution.

 @copyright 2021-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
                      Univ. Bordeaux. All rights reserved.

 @version 6.2.1
 @author Tony Delarue
 @date 2021-10-13

"""
import generation_utils as gu

def genIparmCompletion( iparms ) :
    result = "            COMPREPLY=($(compgen -W \""

    isize  = len("            COMPREPLY=($(compgen -W \"")
    indent = ""
    for group in iparms:
        if group["name"] == "subset_for_old_interface" :
            continue

        for iparm in group["subgroup"]:
            if iparm['access'] != 'IN' :
                continue
            result += indent + iparm["name"] + " \\\n"
            indent  = " " * isize

    result = result[:len(result) -3] + '''" -- $cur))
            ;;
'''
    return result

def genDparmCompletion( dparms ) :
    result = "            COMPREPLY=($(compgen -W \""

    isize  = len("            COMPREPLY=($(compgen -W \"")
    indent = ""
    for dparm in dparms:
        if dparm['access'] != 'IN' :
            continue
        result += indent + dparm["name"] + " \\\n"
        indent  = " " * isize

    result = result[:len(result) -3] + '''" -- $cur))
            ;;
'''
    return result

def genEnumsCompletion( iparms, enums ) :
    result = ""
    isize = len("            COMPREPLY=($(compgen -W \"")
    for group in iparms:
        if group["name"] == "subset_for_old_interface" :
            continue

        for iparm in group["subgroup"]:
            if iparm['access'] != 'IN' :
                continue
            if 'enum' not in iparm :
                continue

            enumname  = iparm['enum']
            result += "        " + iparm["name"] + ''')
            COMPREPLY=($(compgen -W "'''
            values = gu.findEnumFromName( enumname, enums )['values']
            indent = ""
            for value in values :
                name = value["name"]
                result += indent + name.lower() + " \\\n"
                indent  = " " * isize
            result = result[:len(result) -3] + '''" -- $cur))
            ;;
'''
            # -o or --ord
            if iparm["name"] == "iparm_ordering" :
                enumsize = len("pastixorder")
                result += '''        -o|--ord)
            COMPREPLY=($(compgen -W "'''
                indent = ""
                for value in values :
                    name = value["name"]
                    result += indent + name.lower()[enumsize:] + " \\\n"
                    indent  = " " * isize
                result = result[:len(result) -3] + '''" -- $cur))
            ;;
'''
    return result

def genShortCompletion():
    result = '''        -v|--verbose)
            COMPREPLY=($(compgen -W "0 1 2 3" -- $cur))
            ;;
        -f|--fact)
            COMPREPLY=($(compgen -W "0 1 2 3 4" -- $cur))
            ;;
        -s|--sched)
            COMPREPLY=($(compgen -W "0 1 2 3 4" -- $cur))
            ;;
        -c|--check)
            COMPREPLY=($(compgen -W "0 1 2 3 4 5 6" -- $cur))
            ;;
'''
    return result

def genNoCompletion() :
    no_long  = "--rsa --hb --ijv --mm --spm --lap --xlap --graph --threads --gpus --help"
    no_short = "-0 -1 -2 -3 -4 -9 -x -G -t -g -h"

    indent = " " * 8
    result = indent + "# For remaining options, we don't suggest anything\n" + indent

    opts = no_short.split( " " )
    for opt in opts :
        result += opt + "|"
    result += " \\\n" + indent

    linesize = len(no_short)
    size = 0
    opts = no_long.split( " " )
    for opt in opts :
        size   += len(opt) + 1
        result += opt + "|"
        if size > linesize :
            size = 0
            result += " \\\n" + indent

    result += "iparm_*|dparm_*)\n"

    result += indent + '''    COMPREPLY=($(compgen -W "" -- $cur))
            ;;
'''
    return result

def genCompleteCommand():
    result = '''# Add the dynamic completion to the executable
for exec in $(find $BINARY_DIR -maxdepth 1 -executable -type f)
do
    e=$(basename $exec)
    complete -F _pastix_completion $e
done
'''
    return result

def genCompletion( iparms, dparms, enums ) :
    """
    Write the pastix_completion script.

    @in  iparms : The array containing the iparms.
    @in  dparms : The array containing the dparms.
    @in  enums  : The array containing the enums.
    @out  : The complete pastix_completion.bash string.
    """
    isize = len("    local LONG_OPTIONS=(\"") * " "
    long_opts  = '''--rsa --hb --ijv --mm --spm --lap --xlap --graph \\
'''+ isize +'''--threads --gpus --sched --ord --fact --check --iparm --dparm --verbose --help'''
    short_opts = '''-0 -1 -2 -3 -4 -9 -x -G \\
'''+ isize +''' -t -g -s -o -f -c -i -d -v -h'''

    result = '''#
# @file pastix_completion.sh
#
'''+ gu.const_str.replace(" *", "#") +'''
#
#!/usr/bin/env bash

BINARY_DIR=@CMAKE_INSTALL_PREFIX@/examples

_pastix_completion()
{
    local LONG_OPTIONS=("'''+ long_opts +'''")
    local SHORT_OPTIONS=("'''+ short_opts +'''")

    local i cur=${COMP_WORDS[COMP_CWORD]}

    COMPREPLY=($(compgen -W "${LONG_OPTIONS[@]} ${SHORT_OPTIONS[@]}" -- $cur))

    prev=${COMP_WORDS[COMP_CWORD-1]}
    case $prev in
        -i|--iparm)
''' + genIparmCompletion(iparms) + '''
        -d|--dparm)
''' + genDparmCompletion(dparms) + '''
''' + genEnumsCompletion( iparms, enums) + '''
''' + genShortCompletion() + '''
''' + genNoCompletion() +'''
    esac
}

''' + genCompleteCommand()
    return result
