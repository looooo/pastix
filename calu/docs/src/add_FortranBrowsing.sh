cd ../../core_blas
LIST1=`ls -1 core_*.f`
cd ../core_lapack
LIST2=`ls -1 *.f`
cd ../docs/html/htmlbrowsing
for f in `echo $LIST1 $LIST2`
    do
    ROUTINE=`basename $f | sed -e 's/\.f//'`
    sed -i  -e "s/$ROUTINE/<a href=\"$ROUTINE.f.html\">$ROUTINE<\/a>/g" core_*.c.html
done