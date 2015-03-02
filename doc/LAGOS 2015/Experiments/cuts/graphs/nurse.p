set title 'INRC'
set key inside bottom vertical Right
set xlabel 'time (sec)'
set ylabel 'average gap closed (%)'
set xrange [0:300]
set yrange [30:80]
set term postscript eps enhanced "Arial" 24
set output 'nurse.eps'
plot 'nurseSEPEXT.txt' using 1:2 title 'LNPSEP' with linespoints , \
     'nurseSEP.txt' using 1:2 title 'NPSEP' with linespoints , \
     'nurseCGL.txt' using 1:2 title 'CGL' with linespoints
     