#!/bin/bash

ulimit -s 65536

for i in miplib2010-benchmark/*.mps
do
	./cptest $i >> miplib.log
done

for i in INRC/*.mps
do
	./cptest $i >> nurse.log
done

for i in Telebus/*.mps
do
	./cptest $i >> telebus.log
done

for i in Uchoa/*.mps
do
	./cptest $i >> uchoa.log
done