set title 'MIPLIB'
set key inside bottom vertical Right
set xlabel 'time (sec)'
set ylabel 'average gap closed (%)'
set xrange [0:150]
set yrange [30:70]
set term postscript eps enhanced "Arial" 24
set output 'miplib.eps'
plot 'miplibSEPEXT.txt' using 1:2 title 'LNPSEP' with linespoints ls 4, \
     'miplibSEP.txt' using 1:2 title 'NPSEP' with linespoints ls 6, \
     'miplibCGL.txt' using 1:2 title 'CGL' with linespoints ls 8
     