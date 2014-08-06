set title "Lower Bound Improvement: long05"
set xlabel "time (sec)"
set terminal jpeg
set output "long05Perc.jpg"
set ylabel "lb improvement (%)"
plot "long05.gpd" using 3:5 title "npsep" with linespoints,          "long05_cgl.gpd" using 3:5 title "CGL" with linespoints
