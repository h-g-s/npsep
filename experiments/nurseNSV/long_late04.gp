set title "Lower Bound Improvement: long_late04"
set xlabel "time (sec)"
set terminal jpeg
set output "long_late04Perc.jpg"
set ylabel "lb improvement (%)"
plot "long_late04.gpd" using 3:5 title "npsep" with linespoints,          "long_late04_cgl.gpd" using 3:5 title "CGL" with linespoints
