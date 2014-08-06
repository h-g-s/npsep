set title "Lower Bound Improvement: sprint_hidden07"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_hidden07Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_hidden07_j.gpd" using 3:5 title "npsep" with linespoints,          "sprint_hidden07_j_cgl.gpd" using 3:5 title "CGL" with linespoints
