#!/bin/bash

# file with all different improvement times for medium instances
echo -n "" > mediumTimesT.txt

for file in medium*.tlog;
do
   cat $file | cut -d "," -f 2 >> mediumTimesT.txt
done

cat mediumTimesT.txt | sort -n | \
       cut -d "." -f 1 | \
       uniq > mediumTimesT2.txt
rm mediumTimesT.txt

#removing times out of interval
cat mediumTimesT2.txt | grep -v "7[0-9][0-9]" | grep -v "8[0-9][0-9]" > mediumTimes.txt
rm mediumTimesT2.txt

rm cglImprovementmediumAll.txt -f
rm improvementmediumAll.txt -f

echo computing CGL bound improvement:
echo "time"
for cTime in `cat mediumTimes.txt`;
do
   echo -n " $cTime"
   echo -n $cTime >> cglImprovementmediumAll.txt
   for inst in medium*_j_cgl.tlog;
   do
      probName="$(basename $inst _j_cgl.tlog)"
      bbAtTime="$(./bestBoundAtTime.sh $inst $cTime)"
      dist="$(./lbDistance.sh $probName $bbAtTime)"
      echo -n " $dist" >> cglImprovementmediumAll.txt
   done
   
   echo "" >> cglImprovementmediumAll.txt
done

echo ""

echo computing NPSep bound improvement:
echo "time"
for cTime in `cat mediumTimes.txt`;
do
   echo -n " $cTime"
   echo -n $cTime >> improvementmediumAll.txt
   for inst in medium*_j.tlog;
   do
      probName="$(basename $inst _j.tlog)"
      bbAtTime="$(./bestBoundAtTime.sh $inst $cTime)"
      dist="$(./lbDistance.sh $probName $bbAtTime)"
      echo -n " $dist" >> improvementmediumAll.txt
   done
   
   echo "" >> improvementmediumAll.txt
done

