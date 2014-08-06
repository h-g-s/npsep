#!/bin/bash

if [ $# -ne 2 ]
then
   echo parameters: instance_name  time
   exit 0
fi

instName="$1"
tim="$2"
bestObj="0"

for line in `cat $instName`;
do
#   echo $line
   obj="`echo $line | cut -d "," -f 1`"
   currTime="`echo $line | cut -d "," -f 2`"
   stillInTime=$(echo "${currTime}<=${tim}" | bc)
   if [ "$stillInTime" -eq "1" ]
   then
      bestObj="$obj"
   fi
#   echo obj: $obj  time: $currTime in time $stillInTime
done

echo $bestObj

