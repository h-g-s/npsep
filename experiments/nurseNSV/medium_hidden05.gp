set title "Lower Bound Improvement: medium_hidden05"
set xlabel "time (sec)"
set terminal jpeg
set output "medium_hidden05Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium_hidden05.gpd" using 3:5 title "npsep" with linespoints,          "medium_hidden05_cgl.gpd" using 3:5 title "CGL" with linespoints
