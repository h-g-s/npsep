#!/bin/bash

for file in *_log.txt;
do
   pName="`echo $file | sed -e 's/_log.txt//g'`"
   ../logParser/bin/dbg/logParser $file > ${pName}.gpd
done
