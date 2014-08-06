set title "Lower Bound Improvement: sprint_hidden06"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_hidden06Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_hidden06_j.gpd" using 3:5 title "npsep" with linespoints,          "sprint_hidden06_j_cgl.gpd" using 3:5 title "CGL" with linespoints
