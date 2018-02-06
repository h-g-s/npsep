for file in ~/inst/cbcbench/*.mps.gz;
do
    inst=`basename $file .mps.gz`
    echo -n "$inst," >> summary.csv
    ./bin/opt/mrgclq $file 2 1 1
done
for file in ~/inst/bpconf/lp/*.lp;
do
    inst=`basename $file .lp`
    echo -n "$inst,"  >> summary.csv
    ./bin/opt/mrgclq $file 2 1 1
done

