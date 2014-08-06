#!/bin/bash

for file in ../../instances/mips/*.lp;
do
   echo running for $file
   pName=`basename $file .lp`
   echo $pName
   ./writeclqw "${file}" > ${pName}.txt
done

