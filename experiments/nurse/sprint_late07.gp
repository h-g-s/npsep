set title "Lower Bound Improvement: sprint_late07"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_late07Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_late07_j.gpd" using 3:5 title "npsep" with linespoints,          "sprint_late07_j_cgl.gpd" using 3:5 title "CGL" with linespoints
