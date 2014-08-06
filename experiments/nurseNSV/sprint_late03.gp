set title "Lower Bound Improvement: sprint_late03"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_late03Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_late03.gpd" using 3:5 title "npsep" with linespoints,          "sprint_late03_cgl.gpd" using 3:5 title "CGL" with linespoints
