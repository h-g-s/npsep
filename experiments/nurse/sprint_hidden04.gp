set title "Lower Bound Improvement: sprint_hidden04"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_hidden04Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_hidden04_j.gpd" using 3:5 title "npsep" with linespoints,          "sprint_hidden04_j_cgl.gpd" using 3:5 title "CGL" with linespoints
