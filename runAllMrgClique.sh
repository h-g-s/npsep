rm summary.csv
for file in ~/inst/cbcbench/*.mps.gz;
do
    inst=`basename $file .mps.gz`
    echo -n "$inst," >> summary.csv
    ./Release/mrgclq $file 2 1 1
done
for file in ~/inst/bpconf/lp/*.lp;
do
    inst=`basename $file .lp`
    echo -n "$inst,"  >> summary.csv
    ./Release/mrgclq $file 2 1 1
done
mv summary.csv summary-2clq.csv
for file in ~/inst/cbcbench/*.mps.gz;
do
    inst=`basename $file .mps.gz`
    echo -n "$inst," >> summary.csv
    ./Release/mrgclq $file 1 1 1
done
for file in ~/inst/bpconf/lp/*.lp;
do
    inst=`basename $file .lp`
    echo -n "$inst,"  >> summary.csv
    ./Release/mrgclq $file 1 1 1
done
mv summary.csv summary-1clq.csv
for file in ~/inst/cbcbench/*.mps.gz;
do
    inst=`basename $file .mps.gz`
    echo -n "$inst," >> summary.csv
    ./Release/mrgclq $file 3 1 1
done
for file in ~/inst/bpconf/lp/*.lp;
do
    inst=`basename $file .lp`
    echo -n "$inst,"  >> summary.csv
    ./Release/mrgclq $file 3 1 1
done
mv summary.csv summary-3clq.csv

