set title "Lower Bound Improvement: sprint_hidden09"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_hidden09Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_hidden09_j.gpd" using 3:5 title "npsep" with linespoints,          "sprint_hidden09_j_cgl.gpd" using 3:5 title "CGL" with linespoints
