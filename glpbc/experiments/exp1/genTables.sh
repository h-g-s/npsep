#!/bin/bash

timeLimit=`cat runTests.sh | grep -i timeLimit= | cut -d "=" -f 2`
minViol=`cat runTests.sh | grep -i minViol= | cut -d "=" -f 2`
maxDepth=`cat runTests.sh | grep -i maxDepth= | cut -d "=" -f 2`
maxPasses=`cat runTests.sh | grep -i maxPasses= | cut -d "=" -f 2`

echo "- Improvement in Dual Bound -"
echo
echo
echo - timeLimit: $timeLimit
echo - minViol: $minViol
echo - maxDepth: $maxDepth
echo - maxPasses: $maxPasses

echo -n "|| Problem Name "

for useClique in `seq 0 1`;
do
   for useOtherCuts in `seq 0 1`;
   do
      echo -n "| ($useClique  $useOtherCuts) "
   done
done
echo "|"

for file in ~/Dropbox/inst/miplibClique/*.lp;
do
   probName="`basename ${file} .lp`"
   echo -n "| $probName "
   for useClique in `seq 0 1`;
   do
      for useOtherCuts in `seq 0 1`;
      do
         sCCuts="-"
         sOCuts="-"

         if [ $useClique = "1" ]; then
            sCCuts="cCuts"
         fi

         if [ $useOtherCuts = "1" ]; then
            sOCuts="oCuts"
         fi

         logFN="$probName-log-dbound-$sCCuts-$sOCuts-$maxDepth-$minViol-$maxPasses.txt"
         bestBound=`cat $logFN | grep "^$timeLimit" | cut -d " " -f 2`
         sBetter=""
         if [ $useClique = "1" ]; then
            if [ $useOtherCuts = "0" ]; then
               logFNNoCliques="$probName-log-dbound-----$maxDepth-$minViol-$maxPasses.txt"
               bestBoundNoCliques=`cat $logFNNoCliques | grep "^$timeLimit" | cut -d " " -f 2`
               if [ `echo $bestBoundNoCliques \< $bestBound | bc` = "1" ]; then
                  sBetter="*"
               fi
            fi
         fi
         echo -n "|     $bestBound$sBetter "

      done
   done

   echo \|
done
