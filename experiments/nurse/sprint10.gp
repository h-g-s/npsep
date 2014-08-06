set title "Lower Bound Improvement: sprint10"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint10Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint10_j.gpd" using 3:5 title "npsep" with linespoints,          "sprint10_j_cgl.gpd" using 3:5 title "CGL" with linespoints
