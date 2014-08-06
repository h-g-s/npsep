set title "Lower Bound Improvement: medium_hidden02"
set xlabel "time (sec)"
set terminal jpeg
set output "medium_hidden02Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium_hidden02.gpd" using 3:5 title "npsep" with linespoints,          "medium_hidden02_cgl.gpd" using 3:5 title "CGL" with linespoints
