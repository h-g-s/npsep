set title "Lower Bound Improvement: long_hint03"
set xlabel "time (sec)"
set terminal jpeg
set output "long_hint03Perc.jpg"
set ylabel "lb improvement (%)"
plot "long_hint03.gpd" using 3:5 title "npsep" with linespoints,          "long_hint03_cgl.gpd" using 3:5 title "CGL" with linespoints
