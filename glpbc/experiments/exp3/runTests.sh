timeLimit=600
minViol=1050
maxDepth=1
maxPasses=3

for file in ~/Dropbox/inst/miplibClique/*.lp;
do
   echo $file
   echo $timeLimit

   for useClique in `seq 0 1`;
   do
      for useOtherCuts in `seq 0 1`;
      do
         ../glpbc \
            $file $useClique $useOtherCuts $timeLimit $maxDepth $minViol $maxPasses 1
      done
   done
done

