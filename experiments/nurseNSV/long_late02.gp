set title "Lower Bound Improvement: long_late02"
set xlabel "time (sec)"
set terminal jpeg
set output "long_late02Perc.jpg"
set ylabel "lb improvement (%)"
plot "long_late02.gpd" using 3:5 title "npsep" with linespoints,          "long_late02_cgl.gpd" using 3:5 title "CGL" with linespoints
