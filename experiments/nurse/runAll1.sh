#!/bin/bash
for file in /home/haroldo/dev/nurse/instances/*_j.lp;
do
   echo runing for problem: $file
   probName=`basename $file .lp`
   ./cbcr $file -maxPasses=30 | tee ${probName}_log.txt
done

