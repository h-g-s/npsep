set title "Lower Bound Improvement: long03"
set xlabel "time (sec)"
set terminal jpeg
set output "long03Perc.jpg"
set ylabel "lb improvement (%)"
plot "long03_j.gpd" using 3:5 title "npsep" with linespoints,          "long03_j_cgl.gpd" using 3:5 title "CGL" with linespoints
