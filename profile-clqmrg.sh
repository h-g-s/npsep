CPUPROFILE=cpuprofile.log ./bin/dbg/mrgclq ~/inst/cbcbench/neos-506428.mps.gz 2 1 1
pprof --gv ./bin/dbg/mrgclq cpuprofile.log
