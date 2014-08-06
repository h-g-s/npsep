set title "Lower Bound Improvement: sprint_late06"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_late06Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_late06.gpd" using 3:5 title "npsep" with linespoints,          "sprint_late06_cgl.gpd" using 3:5 title "CGL" with linespoints
