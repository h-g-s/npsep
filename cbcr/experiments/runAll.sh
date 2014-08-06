#!/bin/bash
for file in /home/haroldo/dev/npsep/instances/mips/*.lp;
do
   echo runing for problem: $file
   probName=`basename $file .lp`
   ./cbcr $file -maxPasses 5 | tee ${probName}_log.txt
   ./cbcr $file -cgl -maxPasses 5 | tee ${probName}_cgl_log.txt
done

