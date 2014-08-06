set title "Lower Bound Improvement: long01"
set xlabel "time (sec)"
set terminal jpeg
set output "long01Perc.jpg"
set ylabel "lb improvement (%)"
plot "long01_j.gpd" using 3:5 title "npsep" with linespoints,          "long01_j_cgl.gpd" using 3:5 title "CGL" with linespoints
