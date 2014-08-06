set title "Lower Bound Improvement: sprint04"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint04Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint04.gpd" using 3:5 title "npsep" with linespoints,          "sprint04_cgl.gpd" using 3:5 title "CGL" with linespoints
