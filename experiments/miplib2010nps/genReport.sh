#!/bin/bash

for file in *.log;
do
   instName=`basename $file .log`
   lastLine=`cat $file | grep -i "end of root node"`
   timeStr=`echo "$lastLine" | cut -d ":" -f 3 | cut -d " " -f 2`
   echo $instName,$timeStr
done

