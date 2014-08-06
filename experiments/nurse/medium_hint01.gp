set title "Lower Bound Improvement: medium_hint01"
set xlabel "time (sec)"
set terminal jpeg
set output "medium_hint01Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium_hint01_j.gpd" using 3:5 title "npsep" with linespoints,          "medium_hint01_j_cgl.gpd" using 3:5 title "CGL" with linespoints
