#!/bin/bash

for file in ../../instances/mips/*.lp;
do
   pName=`basename $file .lp`
   echo $pName
   ./writeclqw $file > $pName.txt
done

