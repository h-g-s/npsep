#!/bin/bash

if [ "$#" -eq "0" ]
then
   echo usage: instance lb
   exit 0
fi

inst=$1
lb=$2

#echo inst: $inst
#echo lb: $lb

#ub="`cat ub.txt | grep -i '${inst};' | cut -d ';' -f 2`"
ub=$(cat ub.txt | grep -i "${inst};" | cut -d ";" -f 2)

dist=$( echo " scale=4 ; 100*((${ub}-${lb})/${ub})" | bc )

if [ -z "$ub" ]
then
   echo ERROR ERROR, INSTANCE $inst not found
   echo ERROR ERROR, INSTANCE $inst not found
   echo ERROR ERROR, INSTANCE $inst not found
   exit 1
fi

#echo ub: $ub
echo $dist

