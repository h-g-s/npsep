#!/bin/bash
timeLimit=600
minViol=1050
maxDepth=3
maxPasses=3
useClique=1
useOtherCuts=0
problemsDir=~/Dropbox/inst/miplibClique

function run_glpbc 
{
   echo ../glpbc $file $useClique $useOtherCuts $timeLimit $maxDepth $minViol $maxPasses 1
}

function back_to_default_parameters
{
   timeLimit=600
   minViol=1050
   maxDepth=3
   maxPasses=3
   useClique=1
   useOtherCuts=0
}

# start:
for file in $problemsDir/*.lp;
do
   run_glpbc
done

