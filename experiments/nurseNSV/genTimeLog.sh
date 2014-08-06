#!/bin/bash
for file in *_log.txt;
do
   echo runing for problem: $file
   probName=`basename $file .lp`
   tlogFile="`echo $file | sed -e 's/_log.txt/\.tlog/g'`"
   cat $file | grep -i round | sed -e 's/time:/,/g' | cut -d ":" -f 3 | sed -e 's/ //g' > $tlogFile
done

