#!/bin/sh

if [ $# -gt 0 ]
then
    BUILDDIR=$1
fi
BUILDDIR=${BUILDDIR-=build}

SRCDIR_TO_ANALYZE="$BUILDDIR bcsc blend common example graph include kernels order refinement sopalin spm symbol test"

echo $PWD
rm -f filelist.txt

git ls-files | grep "\.[ch]"   >  filelist.txt
git ls-files | grep "\.py"     >> filelist.txt
find $BUILDDIR -name '*\.[ch]' >> filelist.txt
echo "wrappers/python/examples/pypastix/enum.py" >> filelist.txt

# Remove files in kernel/gpus that are C++ and not our own files.
sed -i "/kernels\/gpus\/.*/d" filelist.txt

# Remove all CMakeFiles generated file
sed -i '/CMakeFiles/d' filelist.txt

# Remove files compiled from jdf files
for jdf in `find -name "*\.jdf"`
do
    name=`basename $jdf .jdf`
    echo Remove $name
    sed -i "/$name\.[ch]/d" filelist.txt
done

# Remove installed files
sed -i '/^install.*/d' filelist.txt

# Remove original files used for precision generation
for file in `git grep "@precisions" | awk -F ":" '{ print $1 }'`
do
    sed -i "\:^$file.*:d" filelist.txt
done

# Remove external header files
for file in include/cblas.h include/lapacke.h
do
    sed -i "\:^$file.*:d" filelist.txt
done

# Remove external driver files
for file in spm/drivers/iohb.c spm/drivers/iohb.h spm/drivers/mmio.c spm/drivers/mmio.h
do
    sed -i "\:^$file.*:d" filelist.txt
done

sed '/\.h$/d' filelist.txt > filelist-c.txt

