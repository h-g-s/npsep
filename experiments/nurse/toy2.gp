set title "Lower Bound Improvement: toy2"
set xlabel "time (sec)"
set terminal jpeg
set output "toy2Perc.jpg"
set ylabel "lb improvement (%)"
plot "toy2_j.gpd" using 3:5 title "npsep" with linespoints,          "toy2_j_cgl.gpd" using 3:5 title "CGL" with linespoints
