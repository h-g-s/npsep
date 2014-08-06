set title "Lower Bound Improvement: sprint_hidden02"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_hidden02Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_hidden02.gpd" using 3:5 title "npsep" with linespoints,          "sprint_hidden02_cgl.gpd" using 3:5 title "CGL" with linespoints
