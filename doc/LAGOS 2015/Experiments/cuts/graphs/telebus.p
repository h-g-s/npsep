set title 'Telebus'
set key inside bottom vertical Right
set xlabel 'time (sec)'
set ylabel 'average gap closed (%)'
set xrange [0:300]
set yrange [30:80]
set term postscript eps enhanced "Arial" 24
set output 'telebus.eps'
plot 'telebusSEPEXT.txt' using 1:2 title 'LNPSEP' with linespoints , \
	 'telebusSEP.txt' using 1:2 title 'NPSEP' with linespoints , \
     'telebusCGL.txt' using 1:2 title 'CGL' with linespoints
     