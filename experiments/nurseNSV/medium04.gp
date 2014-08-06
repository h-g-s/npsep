set title "Lower Bound Improvement: medium04"
set xlabel "time (sec)"
set terminal jpeg
set output "medium04Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium04.gpd" using 3:5 title "npsep" with linespoints,          "medium04_cgl.gpd" using 3:5 title "CGL" with linespoints
