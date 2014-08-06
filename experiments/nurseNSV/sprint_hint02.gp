set title "Lower Bound Improvement: sprint_hint02"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_hint02Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_hint02.gpd" using 3:5 title "npsep" with linespoints,          "sprint_hint02_cgl.gpd" using 3:5 title "CGL" with linespoints
