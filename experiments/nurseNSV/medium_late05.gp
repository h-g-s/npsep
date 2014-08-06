set title "Lower Bound Improvement: medium_late05"
set xlabel "time (sec)"
set terminal jpeg
set output "medium_late05Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium_late05.gpd" using 3:5 title "npsep" with linespoints,          "medium_late05_cgl.gpd" using 3:5 title "CGL" with linespoints
