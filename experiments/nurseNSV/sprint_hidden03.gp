set title "Lower Bound Improvement: sprint_hidden03"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_hidden03Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_hidden03.gpd" using 3:5 title "npsep" with linespoints,          "sprint_hidden03_cgl.gpd" using 3:5 title "CGL" with linespoints
