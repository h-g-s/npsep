#!/bin/bash 
n="`cat $1 | grep -i "original cgra" | cut -d ":" -f 3 | cut -d " " -f 2`"
m="`cat $1 | grep -i "original cgra" | cut -d ":" -f 4 | cut -d " " -f 2`"
t="`cat $1 | grep -i "original cgra" | cut -d ":" -f 5 | cut -d " " -f 2`"
p="`basename $1 .txt`"
r="`cat $1 |  grep variables | grep rows | cut -d ")" -f 2 | cut -d " " -f 2`"
nz=`cat $1 | grep -i "variables" | grep -i "nz" | cut -d "s" -f 3 | cut -d " " -f 2`
echo "$p;$n;$m;$r;$nz;$t"
