set title "Lower Bound Improvement: medium_hidden04"
set xlabel "time (sec)"
set terminal jpeg
set output "medium_hidden04Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium_hidden04.gpd" using 3:5 title "npsep" with linespoints,          "medium_hidden04_cgl.gpd" using 3:5 title "CGL" with linespoints
