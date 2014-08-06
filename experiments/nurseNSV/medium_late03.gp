set title "Lower Bound Improvement: medium_late03"
set xlabel "time (sec)"
set terminal jpeg
set output "medium_late03Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium_late03.gpd" using 3:5 title "npsep" with linespoints,          "medium_late03_cgl.gpd" using 3:5 title "CGL" with linespoints
