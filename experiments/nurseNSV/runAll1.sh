#!/bin/bash
for file in /home/haroldo/dev/nurse/instances/*_j.lp;
do
   probName=`basename $file _j.lp`
   echo runing for problem: $probName
   dir=`pwd`
   cd /tmp 
   ~/bin/gmp2lp.sh ~/dev/nurse/modelo/nurseJanNSV.mod ~/dev/nurse/instances/${probName}_j.dat 
   cd $dir
   ./cbcr /tmp/${probName}_j.lp  -maxPasses=30 | tee ${probName}_log.txt
   ./cbcr /tmp/${probName}_j.lp -cgl  -maxPasses=30 | tee ${probName}_cgl_log.txt
done

