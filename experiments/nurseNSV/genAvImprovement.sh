#!/bin/bash

if [ "$#" -eq "0" ]
then
   echo usage: genAvImprovement fileWithAllImprovementsPerTime
   exit 0
fi

dataFile=$1

cat $dataFile |
(
   while read line
   do
#echo line

      sum="0"

      strTime=$(echo $line | cut -d " " -f 1)
      strDists=$(echo $line | cut -d " " -f 2-)

      nEl=0

      for exp in $strDists;
      do
#         echo exp: $exp
         sum=$(echo "scale=5 ; ${sum}+${exp}" | bc )
         nEl=$(($nEl+1))
      done

      av=$(echo "scale=2 ; ${sum}/${nEl}" | bc)
      echo $strTime $av
#echo line sum is $sum $nEl $av
   done
)

