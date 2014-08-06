set title "Lower Bound Improvement: sprint_hidden01"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_hidden01Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_hidden01.gpd" using 3:5 title "npsep" with linespoints,          "sprint_hidden01_cgl.gpd" using 3:5 title "CGL" with linespoints
