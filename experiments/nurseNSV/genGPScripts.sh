#!/bin/bash


for file in *cgl_log.txt;
do
   fName=`basename $file`
   pName=`echo $fName | sed -e 's/_cgl_log.txt//g'`
   echo running for $pName
   logNPS="${pName}.gpd"
   logCGL="${pName}_cgl.gpd"
   gpScr="${pName}.gp"

   echo "set title \"Lower Bound Improvement: $pName\"" > $gpScr
   echo "set xlabel \"time (sec)\"" >> $gpScr
   echo "set terminal jpeg" >> $gpScr
   echo "set output \"${pName}Perc.jpg\"" >> $gpScr
   echo "set ylabel \"lb improvement (%)\"" >> $gpScr
   echo "plot \"${logNPS}\" using 3:5 title \"npsep\" with linespoints, \
         \"${logCGL}\" using 3:5 title \"CGL\" with linespoints" >> $gpScr

   gnuplot $gpScr
done
