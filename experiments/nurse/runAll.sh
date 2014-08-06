#!/bin/bash
for file in /home/haroldo/dev/nurse/instances/*_j.lp;
do
   echo runing for problem: $file
   probName=`basename $file .lp`
   ./cbcr $file -maxPasses=10 | tee ${probName}_log.txt
   ./cbcr $file -cgl -maxPasses=10 | tee ${probName}_cgl_log.txt
done

