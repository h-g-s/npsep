set title "Lower Bound Improvement: sprint_late09"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_late09Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_late09.gpd" using 3:5 title "npsep" with linespoints,          "sprint_late09_cgl.gpd" using 3:5 title "CGL" with linespoints
