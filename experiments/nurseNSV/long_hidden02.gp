set title "Lower Bound Improvement: long_hidden02"
set xlabel "time (sec)"
set terminal jpeg
set output "long_hidden02Perc.jpg"
set ylabel "lb improvement (%)"
plot "long_hidden02.gpd" using 3:5 title "npsep" with linespoints,          "long_hidden02_cgl.gpd" using 3:5 title "CGL" with linespoints
