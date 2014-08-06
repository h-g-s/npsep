set title "Lower Bound Improvement: long_hint02"
set xlabel "time (sec)"
set terminal jpeg
set output "long_hint02Perc.jpg"
set ylabel "lb improvement (%)"
plot "long_hint02_j.gpd" using 3:5 title "npsep" with linespoints,          "long_hint02_j_cgl.gpd" using 3:5 title "CGL" with linespoints
