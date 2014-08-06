set title "Lower Bound Improvement: sprint_hidden08"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_hidden08Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_hidden08.gpd" using 3:5 title "npsep" with linespoints,          "sprint_hidden08_cgl.gpd" using 3:5 title "CGL" with linespoints
