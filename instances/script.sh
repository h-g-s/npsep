#!/bin/bash

ulimit -s 65536

for i in MIPLIB/*.mps.gz
do
	./writeclqw $i >> miplib.log
done

for i in INRC/*.mps.gz
do
	./writeclqw $i >> nurse.log
done
