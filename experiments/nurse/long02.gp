set title "Lower Bound Improvement: long02"
set xlabel "time (sec)"
set terminal jpeg
set output "long02Perc.jpg"
set ylabel "lb improvement (%)"
plot "long02_j.gpd" using 3:5 title "npsep" with linespoints,          "long02_j_cgl.gpd" using 3:5 title "CGL" with linespoints
