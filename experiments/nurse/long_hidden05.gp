set title "Lower Bound Improvement: long_hidden05"
set xlabel "time (sec)"
set terminal jpeg
set output "long_hidden05Perc.jpg"
set ylabel "lb improvement (%)"
plot "long_hidden05_j.gpd" using 3:5 title "npsep" with linespoints,          "long_hidden05_j_cgl.gpd" using 3:5 title "CGL" with linespoints
