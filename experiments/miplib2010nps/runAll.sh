#!/bin/bash

for file in ~/inst/miplib2010/*.mps;
do
   instName=`basename $file .mps`
   echo $file 
   ./cbcr $file -mtm=0 -mgm=0 -maxPasses=9999 -extendC=0 \
        > ${instName}.log
done

