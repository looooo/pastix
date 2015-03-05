#!/bin/bash

rm -f res
touch res

s1=0
s2=1000
ok=0
zero=0
for ((j=0;j<=20;j=j+10));
do
    if [ "$j" -eq "$s1" ]
    then
        j=1000
    fi

    awk '{if ($1 ~ /int/ && $2 ~ /AWK_STOP/) printf "  %s %s = '$j';\n", $1, $2; else print $0}' symbol/symbol_reordering.c > modif
    cat modif > symbol/symbol_reordering.c

    cd build
    make analyze
    cd ..

    echo "\begin{table}[!htbp]" >> res
    echo "\begin{center}" >> res
    echo "\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}" >> res
    echo "\hline" >> res
    echo "Matrix & Nb\_blocks/Cblknbr & Blocks\_size & Nb\_blocks/Cblknbr & Blocks\_size & Time ordering & Time re-ordering & Time with symbolic \\\\" >> res
    echo "\hline" >> res


    for ((i=10;i<=30;i=i+10));
    do

        gmk_m3 $i $i $i grid$i.grf -ggrid$i.xyz
        ./build/example/analyze -G ./grid$i.grf > out.txt

        echo -n "grid $i & " >> res
        echo -n `cat out.txt | grep "&"                                      | awk '/&/ {div = $10 / $2 ; print div " & " $16 " & "}'` >> res
        echo -n `cat out.txt | grep "Time to compute ordering"               | awk '/Time/ {print $5 " & "}'`                  >> res
        echo -n `cat out.txt | grep "TIME TO COMPUTE NEW ORDERING"           | awk '/TIME/ {print $6 " & "}'`                  >> res
        cat out.txt | grep "TIME FOR SYMBOLIC FACTORIZATION AND RE-ORDERING" | awk '/TIME/ {print $7 "\\\\"}'                  >> res

    done

    for i in small oilpan shipsec5 inline audi
    do
        ./build/example/analyze -0 ../../matrix/$i.rsa > out.txt

        echo -n "$i & " >> res
        echo -n `cat out.txt | grep "&"                                      | awk '/&/ {div = $10 / $2 ; print div " & " $16 " & "}'` >> res
        echo -n `cat out.txt | grep "Time to compute ordering"               | awk '/Time/ {print $5 " & "}'`                  >> res
        echo -n `cat out.txt | grep "TIME TO COMPUTE NEW ORDERING"           | awk '/TIME/ {print $6 " & "}'`                  >> res
        cat out.txt | grep "TIME FOR SYMBOLIC FACTORIZATION AND RE-ORDERING" | awk '/TIME/ {print $7 "\\\\"}'                  >> res
    done

    echo "\hline" >> res
    echo "\end{tabular}" >> res
    if [ "$j" -eq "$s2" ]
    then
        echo "\caption{SANS CRITERE DE STOP}" >> res
    else
        echo "\caption{STOP at $j}" >> res
    fi
    echo "\end{center}" >> res
    echo "\end{table}" >> res

    if [ "$ok" -eq "$zero" ]
    then
        echo "\newpage" >> res
        ok=1
    else
        ok=0
    fi


    if [ "$j" -eq "$s2" ]
    then
        j=0
    fi
done

mv res /home/gregoire/PaStiX/RESULTS/tables.tex
cd /home/gregoire/PaStiX/RESULTS
mv scotch_analyse.pdf scotch_analyse_saved.pdf
make
evince scotch_analyse.pdf