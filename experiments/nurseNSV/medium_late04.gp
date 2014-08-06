set title "Lower Bound Improvement: medium_late04"
set xlabel "time (sec)"
set terminal jpeg
set output "medium_late04Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium_late04.gpd" using 3:5 title "npsep" with linespoints,          "medium_late04_cgl.gpd" using 3:5 title "CGL" with linespoints
