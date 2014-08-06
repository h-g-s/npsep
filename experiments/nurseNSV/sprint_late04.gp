set title "Lower Bound Improvement: sprint_late04"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_late04Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_late04.gpd" using 3:5 title "npsep" with linespoints,          "sprint_late04_cgl.gpd" using 3:5 title "CGL" with linespoints
