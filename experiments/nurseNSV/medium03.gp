set title "Lower Bound Improvement: medium03"
set xlabel "time (sec)"
set terminal jpeg
set output "medium03Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium03.gpd" using 3:5 title "npsep" with linespoints,          "medium03_cgl.gpd" using 3:5 title "CGL" with linespoints
