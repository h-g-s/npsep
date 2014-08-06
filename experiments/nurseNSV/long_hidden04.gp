set title "Lower Bound Improvement: long_hidden04"
set xlabel "time (sec)"
set terminal jpeg
set output "long_hidden04Perc.jpg"
set ylabel "lb improvement (%)"
plot "long_hidden04.gpd" using 3:5 title "npsep" with linespoints,          "long_hidden04_cgl.gpd" using 3:5 title "CGL" with linespoints
