set title "Lower Bound Improvement: sprint02"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint02Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint02_j.gpd" using 3:5 title "npsep" with linespoints,          "sprint02_j_cgl.gpd" using 3:5 title "CGL" with linespoints
