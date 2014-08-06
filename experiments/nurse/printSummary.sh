#!/bin/bash

# columns:
#             round 1
# problem   cuts e.time bound 


for file in ~/Dropbox/inst/miplibClique/*.lp;
do
   probName=`basename $file .lp`
   probNameP=`echo -n ${probName} | sed s/_/-/g`
   
   strMin=`cat $file | grep -i minimize | cut -d " " -f 1`

   echo -n "${probNameP} & "
   if [ -n "${strMin}" ]; then
      echo -n "Min & "
   else
      echo -n "Max & "
   fi
   log=${probName}_cgl_log.txt

   lastRound=`cat ${log} | grep -i round | tac | head -1 | cut -d " " -f 2`

   i=1
   cuts=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 6`
   time=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 9`
   dual=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 14`
   relaxBound=`cat ${log} | grep -i "Initial dual bound"  | cut -d " " -f 7`
   echo -n "${cuts} & "
   echo -n "${time} & "
   echo -n "${dual}  &"

   i=$lastRound
   echo -n "${lastRound} &"
   cuts=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 6`
   time=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 9`
   dual=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 14`
   relaxBound=`cat ${log} | grep -i "Initial dual bound"  | cut -d " " -f 7`
   echo -n "${cuts} & "
   echo -n "${time} & "
   echo -n "${dual}  &"

   log=${probName}_log.txt
   lastRound=`cat ${log} | grep -i round | tac | head -1 | cut -d " " -f 2`
   i=1

   cuts=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 6`
   time=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 9`
   dual=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 14`
   relaxBound=`cat ${log} | grep -i "Initial dual bound"  | cut -d " " -f 7`
   echo -n "${cuts} & "
   echo -n "${time} & "
   echo -n "${dual}  &"

   i=$lastRound
   echo -n "${lastRound} &"
   cuts=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 6`
   time=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 9`
   dual=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 14`
   relaxBound=`cat ${log} | grep -i "Initial dual bound"  | cut -d " " -f 7`
   echo -n "${cuts} & "
   echo -n "${time} & "
   echo -n "${dual} "
   echo  '\\\\'

done

