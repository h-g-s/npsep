#!/bin/bash

# columns:
#             round 1
# problem   cuts e.time bound 

objSense="EMPTY"

# info for each technique
# considering  FR - first round
# and LR - lastRound

cutsFR=""
timeFR=""
dualFR=""

dualLRCGL=""
dualLRNPSep=""

npSepBetter="0"
cglBetter="0"
objSense=""
objS=""

getObjSense()
{
   strMin=`cat $file | grep -i minimize | cut -d " " -f 1`
   
   if [ -n "${strMin}" ]; then
      objSense="Min"
      objS="-"
   else
      objSense="Max"
      objS="+"
   fi
}

getLogInfo() {
   getObjSense
   lastRound=`cat ${log} | grep -i round | tac | head -1 | cut -d " " -f 2`

# First round:
   i=1
   cutsFR=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 6`
   timeFR=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 9`
   dualFR=`cat ${log} | grep -i "cuts of pass ${i}" | cut -d " " -f 16`

   i=$lastRound
   cutsLR=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 6`
   timeLR=`cat ${log} | grep -i "round ${i}" | cut -d " " -f 9`
   dualLR=`cat $log | grep -i "end of root" | cut -d " " -f 8`
#   relaxBound=`cat ${log} | grep -i "Initial dual bound"  | cut -d " " -f 7`
}

checkObjNPSepBetter()
{
   npSepBetter="0"

   if test "$objSense"=="Min" 
   then
      if [ "`echo "${dualLRNPSep} > ${dualLRCGL}" | bc`" -eq "1" ]; then
         npSepBetter="1"
      else
         if [ "`echo "${dualLRNPSep} < ${dualLRCGL}" | bc`" -eq "1" ]; then
            npSepBetter="1"
         fi
      fi
   fi
}

checkObjCGLBetter()
{
   cglBetter="0"

   if test "$objSense"=="Min" 
   then
      if [ "`echo "${dualLRCGL} > ${dualLRNPSep}" | bc`" -eq "1" ]; then
         cglBetter="1"
      else
         if [ "`echo "${dualLRCGL} < ${dualLRNPSep}" | bc`" -eq "1" ]; then
            cglBetter="1"
         fi
      fi
   fi
}


for file in /home/haroldo/dev/nurse/instances/*_j.lp;
do
   probName=`basename $file .lp`
   probNameP=`echo -n ${probName} | sed s/_/-/g`

   log=${probName}_cgl_log.txt
   getLogInfo 

   dualLRCGL=${dualLR}
   
   echo -n "${probNameP} & "
   echo -n '\'
   echo -n "multicolumn{1}{c|}{"
   echo -n " ${objS} } & "

# just compute dual of other method
   log=${probName}_log.txt
   getLogInfo
   dualLRNPSep=${dualLR}   
   log=${probName}_cgl_log.txt
   getLogInfo

   echo -n "${cutsFR} & ${timeFR} & ${dualFR} & "
   if [ ${cglBetter} -eq "1" ] ; then
      echo -n "${lastRound} & ${cutsLR} & ${timeLR} & "
      echo -n "\\"
      echo -n "textbf{${dualLR}} & "
   else
      echo -n "${lastRound} & ${cutsLR} & ${timeLR} & ${dualLR} & "
   fi
      
   log=${probName}_log.txt
   getLogInfo
   echo -n "${cutsFR} & ${timeFR} & ${dualFR} & "

   checkObjNPSepBetter
   if [ ${npSepBetter} -eq "1" ] ; then
      echo -n "${lastRound} & ${cutsLR} & ${timeLR} & "
      echo -n "\\"
      echo -n "textbf{${dualLR}}"
   else
      echo -n "${lastRound} & ${cutsLR} & ${timeLR} & ${dualLR}"
   fi
   echo ' \\'
   

done


