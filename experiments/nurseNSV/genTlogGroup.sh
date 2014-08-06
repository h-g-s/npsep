#!/bin/bash

# file with all different improvement times for sprint instances
echo -n "" > sprintTimesT.txt

for file in sprint*.tlog;
do
   cat $file | cut -d "," -f 2 >> sprintTimesT.txt
done

cat sprintTimesT.txt | sort -n | \
       cut -d "." -f 1 | \
       uniq > sprintTimesT2.txt
rm sprintTimesT.txt

#removing times out of interval
cat sprintTimesT2.txt | grep -v "7[0-9][0-9]" | grep -v "3[1-9][0-9]" > sprintTimes.txt
rm sprintTimesT2.txt

rm cglImprovementSprintAll.txt -f
rm improvementSprintAll.txt -f

echo computing CGL bound improvement:
echo "time"
for cTime in `cat sprintTimes.txt`;
do
   echo -n " $cTime"
   echo -n $cTime >> cglImprovementSprintAll.txt
   for inst in sprint*_j_cgl.tlog;
   do
      probName="$(basename $inst _j_cgl.tlog)"
      bbAtTime="$(./bestBoundAtTime.sh $inst $cTime)"
      dist="$(./lbDistance.sh $probName $bbAtTime)"
      echo -n " $dist" >> cglImprovementSprintAll.txt
   done
   
   echo "" >> cglImprovementSprintAll.txt
done

echo ""

echo computing NPSep bound improvement:
echo "time"
for cTime in `cat sprintTimes.txt`;
do
   echo -n " $cTime"
   echo -n $cTime >> improvementSprintAll.txt
   for inst in sprint*_j.tlog;
   do
      probName="$(basename $inst _j.tlog)"
      bbAtTime="$(./bestBoundAtTime.sh $inst $cTime)"
      dist="$(./lbDistance.sh $probName $bbAtTime)"
      echo -n " $dist" >> improvementSprintAll.txt
   done
   
   echo "" >> improvementSprintAll.txt
done

