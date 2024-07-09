#
#  @file validate_expand.sh
#
#  @copyright 2016-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.3.2
#  @author Mathieu Faverge
#  @date 2023-07-21
#
# This script compares the output files generated by the DOF expansion
# testings
#
#!/bin/sh

fltname=$1

rm -f *.sort
for i in `ls -1 expand_*.dat`;
do
    sort $i > ${i}.sort
done

for i in 0
do
    for mtxtype in General Symmetric Hermitian
    do
        basefile=expand_${i}_CSC_0_${mtxtype}_${fltname}_sparse_cp.dat.sort

        if [ -f $basefile ]
        then
            echo "-- Dof=$i, $mtxtype -- "

            for baseval in 0 1
            do
                for fmttype in CSC CSR IJV
                do
                    for storage in sparse dense
                    do
                        for comp in cp ucp
                        do
                            if [ -f expand_${i}_${fmttype}_${baseval}_${mtxtype}_${fltname}_${storage}_${comp}.dat.sort ]
                            then
                                echo -n "---- CSC Sparse CP VS $fmttype $storage $comp $baseval: "
                                diff $basefile expand_${i}_${fmttype}_${baseval}_${mtxtype}_${fltname}_${storage}_${comp}.dat.sort | wc -l
                            fi
                        done
                    done
                done
            done
        fi
    done
done

for i in 1
do
    for mtxtype in General Symmetric Hermitian
    do
        for fmttype in CSC CSR IJV
        do
            for baseval in 0 1
            do
                basefile=expand_${i}_${fmttype}_${baseval}_${mtxtype}_${fltname}_sparse_cp.dat.sort

                if [ -f $basefile ]
                then
                    echo "-- Dof=$i, $mtxtype -- "
                    for storage in sparse dense
                    do
                        for comp in cp ucp
                        do
                            if [ -f expand_${i}_${fmttype}_${baseval}_${mtxtype}_${fltname}_${storage}_${comp}.dat.sort ]
                            then
                               echo -n "---- $fmtype Sparse CP VS $fmttype $storage $comp $baseval: "
                               diff $basefile expand_${i}_${fmttype}_${baseval}_${mtxtype}_${fltname}_${storage}_${comp}.dat.sort | wc -l
                            fi
                        done
                    done
                fi
            done
        done
    done
done

