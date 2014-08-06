set title "Lower Bound Improvement: sprint_hint01"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_hint01Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_hint01.gpd" using 3:5 title "npsep" with linespoints,          "sprint_hint01_cgl.gpd" using 3:5 title "CGL" with linespoints
