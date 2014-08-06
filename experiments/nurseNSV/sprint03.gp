set title "Lower Bound Improvement: sprint03"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint03Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint03.gpd" using 3:5 title "npsep" with linespoints,          "sprint03_cgl.gpd" using 3:5 title "CGL" with linespoints
