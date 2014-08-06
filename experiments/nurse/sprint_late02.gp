set title "Lower Bound Improvement: sprint_late02"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_late02Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_late02_j.gpd" using 3:5 title "npsep" with linespoints,          "sprint_late02_j_cgl.gpd" using 3:5 title "CGL" with linespoints
