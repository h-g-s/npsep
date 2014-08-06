#!/bin/bash

for file in ../../instances/miplib/*.mps.gz;
do
   iName=`basename $file .mps.gz`
   valgrind ./writeclqw $file > ${iName}-log.txt 2> ${iName}-errors.txt
done


