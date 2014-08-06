#!/bin/bash

for file in ../../../instances/miplib/*.mps.gz;
do
   instName=`basename $file .mps.gz`
   echo $file 
   ./cbcr $file -mtm=0 -mgm=0 -maxPasses=9999 -extendC=0 \
        > ${instName}.log
done

