set title 'INRC'
set key inside bottom vertical Right
set xlabel 'time (sec)'
set ylabel 'average gap closed (%)'
set xrange [0:150]
set yrange [30:85]
set term postscript eps enhanced "Arial" 24
set output 'nurse.eps'
plot 'nurseSEPEXT.txt' title 'LNPSEP' with linespoints ls 4, \
     'nurseSEP.txt' title 'NPSEP' with linespoints ls 6, \
     'nurseCGL.txt' title 'CGL' with linespoints ls 8
     