set title "Lower Bound Improvement: medium_hint03"
set xlabel "time (sec)"
set terminal jpeg
set output "medium_hint03Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium_hint03.gpd" using 3:5 title "npsep" with linespoints,          "medium_hint03_cgl.gpd" using 3:5 title "CGL" with linespoints
