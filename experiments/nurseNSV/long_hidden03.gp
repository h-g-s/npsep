set title "Lower Bound Improvement: long_hidden03"
set xlabel "time (sec)"
set terminal jpeg
set output "long_hidden03Perc.jpg"
set ylabel "lb improvement (%)"
plot "long_hidden03.gpd" using 3:5 title "npsep" with linespoints,          "long_hidden03_cgl.gpd" using 3:5 title "CGL" with linespoints
