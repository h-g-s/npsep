#!/bin/bash

# file with all different improvement times for long instances
echo -n "" > longTimesT.txt

for file in long*.tlog;
do
   cat $file | cut -d "," -f 2 >> longTimesT.txt
done

cat longTimesT.txt | sort -n | \
       cut -d "." -f 1 | \
       uniq > longTimesT2.txt
rm longTimesT.txt


#removing times out of interval
cat longTimesT2.txt | grep -v "7[0-9][0-9]" | grep -v "6[0-9][1-9]" > longTimes.txt
rm longTimesT2.txt

rm cglImprovementlongAll.txt -f
rm improvementlongAll.txt -f

echo computing CGL bound improvement:
echo "time"
for cTime in `cat longTimes.txt`;
do
   echo -n " $cTime"
   echo -n $cTime >> cglImprovementlongAll.txt
   for inst in long*_j_cgl.tlog;
   do
      probName="$(basename $inst _j_cgl.tlog)"
      bbAtTime="$(./bestBoundAtTime.sh $inst $cTime)"
      dist="$(./lbDistance.sh $probName $bbAtTime)"
      echo -n " $dist" >> cglImprovementlongAll.txt
   done
   
   echo "" >> cglImprovementlongAll.txt
done

echo ""

echo computing NPSep bound improvement:
echo "time"
for cTime in `cat longTimes.txt`;
do
   echo -n " $cTime"
   echo -n $cTime >> improvementlongAll.txt
   for inst in long*_j.tlog;
   do
      probName="$(basename $inst _j.tlog)"
      bbAtTime="$(./bestBoundAtTime.sh $inst $cTime)"
      dist="$(./lbDistance.sh $probName $bbAtTime)"
      echo -n " $dist" >> improvementlongAll.txt
   done
   
   echo "" >> improvementlongAll.txt
done

