set title "Lower Bound Improvement: medium05"
set xlabel "time (sec)"
set terminal jpeg
set output "medium05Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium05.gpd" using 3:5 title "npsep" with linespoints,          "medium05_cgl.gpd" using 3:5 title "CGL" with linespoints
