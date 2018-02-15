CPUPROFILE=cpuprofile.log ./bin/dbg/mrgclq ~/inst/cbcbench/ns1758913.mps.gz 2 1 1
pprof --gv ./bin/dbg/mrgclq cpuprofile.log
