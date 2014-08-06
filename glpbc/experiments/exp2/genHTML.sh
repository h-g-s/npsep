#!/bin/sh
sh genTables.sh > results.txt
txt2tags -t html results.txt 
