set title "Lower Bound Improvement: long_late03"
set xlabel "time (sec)"
set terminal jpeg
set output "long_late03Perc.jpg"
set ylabel "lb improvement (%)"
plot "long_late03.gpd" using 3:5 title "npsep" with linespoints,          "long_late03_cgl.gpd" using 3:5 title "CGL" with linespoints
