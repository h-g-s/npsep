set title "Lower Bound Improvement: medium_hint02"
set xlabel "time (sec)"
set terminal jpeg
set output "medium_hint02Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium_hint02.gpd" using 3:5 title "npsep" with linespoints,          "medium_hint02_cgl.gpd" using 3:5 title "CGL" with linespoints
