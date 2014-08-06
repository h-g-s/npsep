set title "Lower Bound Improvement: long04"
set xlabel "time (sec)"
set terminal jpeg
set output "long04Perc.jpg"
set ylabel "lb improvement (%)"
plot "long04.gpd" using 3:5 title "npsep" with linespoints,          "long04_cgl.gpd" using 3:5 title "CGL" with linespoints
