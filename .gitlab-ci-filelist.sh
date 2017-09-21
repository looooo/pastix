#!/bin/sh
set -x

SRCDIR_TO_ANALYZE="build bcsc blend common example graph include kernels order refinement sopalin spm symbol test"

echo $PWD
rm -f filelist.txt
for dir in ${SRCDIR_TO_ANALYZE}
do
    find $dir -name '*\.[ch]' >> filelist.txt
done

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
sed -i '/build\/install.*/d' filelist.txt

# Remove original files used for precision generation
for file in `git grep "@precisions" | awk -F ":" '{ print $1 }'`
do
    sed -i "\:^$file.*:d" filelist.txt
done

# Remove external header files
for file in include/cblas.h include/lapacke.h include/lapacke_mangling.h
do
    sed -i "\:^$file.*:d" filelist.txt
done

# Remove external driver files
for file in spm/drivers/iohb.c
            spm/drivers/iohb.h
            spm/drivers/mmio.c
            spm/drivers/mmio.h
do
    sed -i "\:^$file.*:d" filelist.txt
done
