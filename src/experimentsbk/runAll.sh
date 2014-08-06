#!/bin/bash
for file in ../../instances/wClique/*.clqw;
do
    echo
    echo running for problem `basename $file .clqw`
    ./BronKerbosch $file
done
