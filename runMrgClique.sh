#!/bin/bash
for nExt in 3 5;
do
    for bkIt in 512 1024 4096 8192;
    do
        rm summary.csv
        for file in ~/inst/cbcbench/*.mps.gz;
        do
            inst=`basename $file .mps.gz`
            echo -n "$inst," >> summary.csv
            ./Release/mrgclq $file $nExt 1 1 $bkIt
        done
        for file in ~/inst/bpconf/lp/*.lp;
        do
            inst=`basename $file .lp`
            echo -n "$inst,"  >> summary.csv
            ./Release/mrgclq $file $nExt 1 1 $bkIt
        done
        mv summary.csv summary-$nExt-$bkIt.csv
    done
done
