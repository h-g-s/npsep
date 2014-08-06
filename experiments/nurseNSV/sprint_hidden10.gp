set title "Lower Bound Improvement: sprint_hidden10"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_hidden10Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_hidden10.gpd" using 3:5 title "npsep" with linespoints,          "sprint_hidden10_cgl.gpd" using 3:5 title "CGL" with linespoints
