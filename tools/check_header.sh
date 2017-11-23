#
#  @file check_header.sh
#
#  @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.0
#  @author Mathieu Faverge
#  @date 2017-06-24
#
# This script check that basic informations is present and correct in
# headers of source files.
#
#!/bin/sh

#
# Check that each header file contains
#
# #ifndef _filename_
# #define _filename_
# ...
# #endif /* _filename_ */
#
files=`git ls-files | grep '\.h' | grep -v "common/sys/atomic-" | grep -v cblas.h | grep -v lapacke.h`

header=1

print_header()
{
    if [ $header -ne 0 ]
    then
        echo "------ $1 --------"
        header=0
    fi
}

check_header_file()
{
    filename=$1
    basename=`basename $filename .in`

    if [ "$basename" != "CMakeLists.txt" ]
    then
        toto=`grep " @file $basename" $filename`
        if [ $? -ne 0 ]
        then
            toto=`grep " @file .*/$basename" $filename`
        fi

        if [ $? -ne 0 ]
        then
            print_header $filename
            echo -n "@file line missing or incorrect:"; grep "@file" $filename; echo ""
        fi
    fi
}

check_header_copyright()
{
    filename=$1
    basename=`basename $filename`

    toto=`grep -E " @copyright [0-9]{4}-2017 Bordeaux INP" $filename`
    if [ $? -ne 0 ]
    then
        toto=`grep -E " @copyright 20[0-9]{2}      Bordeaux INP" $filename`
    fi

    if [ $? -ne 0 ]
    then
        print_header $filename
        echo -n "@copyright line missing or incorrect:"; grep "@copyright" $filename; echo "";
    fi
}

check_header_version()
{
    filename=$1
    basename=`basename $filename`

    toto=`grep -E " @version [0-9]\.[0-9]\.[0-9]" $filename`
    if [ $? -ne 0 ]
    then
        print_header $filename
        echo -n "@version line missing or incorrect:"; grep "@version" $filename; echo "";
    fi
}

check_header_author()
{
    filename=$1
    basename=`basename $filename`

    toto=`grep -E " @author " $filename`
    if [ $? -ne 0 ]
    then
        print_header $filename
        echo "@author line missing";
    fi
}

check_header_date()
{
    filename=$1
    basename=`basename $filename`

    toto=`grep -E " @date [0-9]{4}-[01][0-9]-[0-3][0-9]" $filename`
    if [ $? -ne 0 ]
    then
        print_header $filename
        echo -n "@date line missing or incorrect"; grep "@version" $filename; echo "";
    fi
}

#
# Check that the given source file contains
#
# @file filename
# @copyright
# @version
# @date
#
check_header()
{
    header=1
    check_header_file $1
    check_header_copyright $1
    check_header_version $1
    check_header_author $1
    check_header_date $1
}



for f in $files
do
    filename=`basename $f .h | sed 's/\.h\.in//' | awk '{print tolower($0)}'`
    macro="_${filename}_h_"
    err=0

    toto=`grep "#ifndef .*$macro" $f`
    ret=$?
    err=$((err + ret))

    toto=`grep "#define .*$macro" $f`
    ret=$?
    err=$((err + ret))

    toto=`grep "#endif /\* .*$macro \*/" $f`
    ret=$?
    err=$((err + ret))

    case $err in
    0)
        #echo $f "OK"
        ;;
    *)
        echo "-------- $f / $filename ------------"
        grep "#ifndef" $f
        grep "#define" $f
        grep "#endif"  $f
        ;;
    esac
done


#
# Check headers
#
files=`git ls-files | grep -v "^\." | grep -v ".*\.md" | grep -v LICENSE | grep -v ".*\.cmake" | grep -v "common/sys/atomic-" | grep -v docs/doxygen | grep -v CTest | grep -v cblas.h | grep -v lapacke.h | grep -v kernels/gpus/kepler  | grep -v kernels/gpus/fermi | grep -v test/matrix`

for f in $files
do
    if [ -d $f ]
    then
        continue;
    fi

    check_header $f
done
