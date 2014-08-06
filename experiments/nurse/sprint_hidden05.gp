set title "Lower Bound Improvement: sprint_hidden05"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_hidden05Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_hidden05_j.gpd" using 3:5 title "npsep" with linespoints,          "sprint_hidden05_j_cgl.gpd" using 3:5 title "CGL" with linespoints
