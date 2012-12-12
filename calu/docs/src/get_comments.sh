#!/bin/bash
INCDIR=../../include
SRCDIR=../../src
rm -rf tmp

    for f in `ls $SRCDIR/auxiliary.c $SRCDIR/context.c $SRCDIR/control.c $SRCDIR/descriptor.c $SRCDIR/tile.c $SRCDIR/workspace*.c $SRCDIR/workspace.c $SRCDIR/[czds][gptu]*.c $SRCDIR/zc*.c $SRCDIR/ds*.c`
    do
    #echo "======== $f ========"
    # GET NAME OF THE FILE
    ROUTINE=`basename $f | sed -e 's/\.c//' -e 's/tile/Tile/'`
    FILENAME=`basename $f`
    ROUTINE=`echo $ROUTINE | tr [[:lower:]] [[:upper:]]`
    AUX=`head $f | grep -c -m 1 "auxiliary"`
    # COUNT ROUTINE IN FILE
    COUNT=`grep  -c "P /// U /// R /// P /// O /// S /// E" $f`
    ALLLINES=`wc -l $f | awk '{print $1}'`
    LINES=$ALLLINES
    cat $f > file
    CUT=0

    for ((  i = 1 ;  i <= $COUNT;  i++  ))
    do
    PLASMAROUTINE=`grep -m 1 "^// PLASMA_" file | awk '{print $2}'`
    echo "== $PLASMAROUTINE"
    echo "<a href=\"http://icl.cs.utk.edu/projectsfiles/plasma/html/htmlbrowsing/$FILENAME.html\">$PLASMAROUTINE</a>" >> tmp

    # PURPOSE SECTION
    START_PURPOSE=`grep  -m 1 -n "P /// U /// R /// P /// O /// S /// E" file | sed -e 's/:/ /' | awk '{print $1}'`
    START_PURPOSE=$((START_PURPOSE+1))
    END_PURPOSE=`grep  -m 1 -n "A /// R /// G" file | sed -e 's/:/ /' | awk '{print $1}'`
    END_PURPOSE=$((END_PURPOSE-1))
    CUT_PURPOSE=$((END_PURPOSE-START_PURPOSE+1))
    if [ $CUT_PURPOSE -gt 1 ]
        then
            echo ".Purpose"
            head -n $END_PURPOSE file | tail -n $CUT_PURPOSE |  sed -e 's/\/\/ //g' | sed -e 's/\/\///g' -e "s/ $ROUTINE/*$PLASMAROUTINE*/" -e "s/$PLASMAROUTINE/*$PLASMAROUTINE*/"
    fi

    # INTERFACE SECTION
    echo ""
    echo ".C Bindings"
    echo [source,c]
    echo "----"
    if [[ $AUX -eq 0 ]]
        then
        MPREC=`echo $PLASMAROUTINE | sed -e 's/PLASMA_//' | cut -c 1,2 | tr [[:upper:]] [[:lower:]]`
            if [[ "$MPREC" == "ds" || "$MPREC" == "zc" ]]
            then
                SUF="_m"
            else
                PREC=`echo $PLASMAROUTINE | sed -e 's/PLASMA_//'  | cut -c 1 | tr [[:upper:]] [[:lower:]]`
                SUF="_$PREC"
            fi
            else
                echo "$FILENAME" | grep -i -q 'workspace_'
            if [ $? -eq 0 ] ; then
                # IT IS AN AUXILIARY ALLOC WORKSPACE ROUTINE
                SUF=`echo "$FILENAME" | sed -e 's/workspace//' -e 's/\.c//'`
            else
                # IT IS AN AUXILIARY ROUTINE, SO NO PRECISION AND NO SUFFIX IN INCLUDE AND BINDINGS
                SUF=''
            fi
     fi
    grep -m 1 -i "$PLASMAROUTINE(" $INCDIR/plasma$SUF.h | sed -e 's/;//'
    echo "----"
    echo ""
    echo ".Fortran Bindings"
    echo [source,c]
    echo "----"
    grep  -i "$PLASMAROUTINE(" fortran${SUF}biding.f | sed -e 's/;//'
    echo "----"
    echo ""

    # ARGUMENT SECTION
    START_ARGS=`grep  -m 1 -n "A /// R /// G" file | sed -e 's/:/ /' | awk '{print $1}'`
    START_ARGS=$((START_ARGS+1))
    if [ `grep  -m 1 -c "R /// E /// T ///" file` -gt 1 ]
        then
            END_ARGS=`grep  -m 1 -n "C /// O /// D /// E" file | sed -e 's/:/ /' | awk '{print $1}'`
        else
            END_ARGS=`grep  -m 1 -n "R /// E /// T ///" file | sed -e 's/:/ /' | awk '{print $1}'`
    fi
    END_ARGS=$((END_ARGS-1))
    CUT_ARGS=$((END_ARGS-START_ARGS+1))
    if [ $CUT_ARGS -gt 1 ]
        then
            echo ".Arguments"
            echo [source,c]
            echo "----"
            head -n $END_ARGS file | tail -n $CUT_ARGS | sed -e 's/\/\/ //g' | sed -e 's/\/\///g'
            echo "----"
    fi

    # RETURN VALUE SECTION
    START_RETVAL=`grep  -m 1 -n "R /// E /// T" file | sed -e 's/:/ /' | awk '{print $1}'`
    START_RETVAL=$((START_RETVAL+1))
    END_RETVAL=`grep  -m 1 -n "C /// O /// D /// E" file | sed -e 's/:/ /' | awk '{print $1}'`
    END_RETVAL=$((END_RETVAL-2))
    CUT_RETVAL=$((END_RETVAL-START_RETVAL+1))
    if [ $CUT_ARGS -gt 1 ]
        then
            echo ".Return Value:"
            echo [source,c]
            echo "----"
            head -n $END_RETVAL file | tail -n $CUT_RETVAL | sed -e 's/\/\///g'
            echo "----"
            echo ""
    fi

   # ONLINE BROWSING SECTION
    echo ""
    echo ".Online Browsing"
    echo "Dive into http://icl.cs.utk.edu/projectsfiles/plasma/html/htmlbrowsing/$FILENAME.html[$PLASMAROUTINE]"
    echo ""

    CUT=`grep  -m 1 -n "C /// O /// D /// E" file | sed -e 's/:/ /' | awk '{print $1}'`
    ALLLINES=`wc -l file | awk '{print $1}'` 
    LINES=`expr $ALLLINES - $CUT`
    tail -n  $LINES file > tmp
    mv tmp file
    done
    done
    rm -f file

