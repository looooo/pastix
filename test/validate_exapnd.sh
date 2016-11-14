#!/bin/sh

fltname=$1

rm -f *.sort
for i in `ls -1 *.dat`;
do
    sort $i > ${i}.sort
done

for i in 0 1
do
    for mtxtype in General Symmetric Hermitian
    do
        basefile=${i}_CSC_0_${mtxtype}_${fltname}_sparse_cp.dat.sort

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

                            echo -n "---- CSC Sparse CP VS $fmttype $storage $comp $baseval: "
                            diff $basefile ${i}_${fmttype}_${baseval}_${mtxtype}_${fltname}_${storage}_${comp}.dat.sort | wc -l
                        done
                    done
                done
            done
        fi
    done
done

