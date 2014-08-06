#!/bin/bash

for file in *cgl_log.txt;
do
   pName="`echo $file | sed -e 's/_cgl_log.txt//g'`"
   logCgl="${pName}_cgl.tlog"
   logNPS="${pName}.tlog"
   cglBound="`tac $logCgl | head -1  | cut -d "," -f 1`"
   npsBound="`tac $logNPS | head -1  | cut -d "," -f 1`"
   echo "$pName;$cglBound;$npsBound"
done


