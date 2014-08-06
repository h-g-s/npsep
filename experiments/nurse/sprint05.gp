set title "Lower Bound Improvement: sprint05"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint05Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint05_j.gpd" using 3:5 title "npsep" with linespoints,          "sprint05_j_cgl.gpd" using 3:5 title "CGL" with linespoints
