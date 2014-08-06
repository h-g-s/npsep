set title "Lower Bound Improvement: long_hint01"
set xlabel "time (sec)"
set terminal jpeg
set output "long_hint01Perc.jpg"
set ylabel "lb improvement (%)"
plot "long_hint01.gpd" using 3:5 title "npsep" with linespoints,          "long_hint01_cgl.gpd" using 3:5 title "CGL" with linespoints
