#!/bin/bash

for file in *cgl_log.txt;
do
   pName="`echo $file | sed -e 's/_j_cgl_log.txt//g'`"
   logCgl="${pName}_j_cgl.tlog"
   logNPS="${pName}_j.tlog"
   cglBound="`tac $logCgl | head -1  | cut -d "," -f 1`"
   npsBound="`tac $logNPS | head -1  | cut -d "," -f 1`"
   echo "$pName;$cglBound;$npsBound"
done


