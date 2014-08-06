set title "Lower Bound Improvement: sprint07"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint07Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint07.gpd" using 3:5 title "npsep" with linespoints,          "sprint07_cgl.gpd" using 3:5 title "CGL" with linespoints
