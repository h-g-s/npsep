set title "Lower Bound Improvement: sprint_hint03"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint_hint03Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint_hint03_j.gpd" using 3:5 title "npsep" with linespoints,          "sprint_hint03_j_cgl.gpd" using 3:5 title "CGL" with linespoints
