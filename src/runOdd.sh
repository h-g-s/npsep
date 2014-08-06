#!/bin/bash
for file in ~/dev/nurse/instances/*_j.lp;
do
   pName="`basename $file _j.lp`"
   echo running for $pName
   ./cbcr $file > logOdd${pName}.txt
done

