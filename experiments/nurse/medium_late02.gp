set title "Lower Bound Improvement: medium_late02"
set xlabel "time (sec)"
set terminal jpeg
set output "medium_late02Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium_late02_j.gpd" using 3:5 title "npsep" with linespoints,          "medium_late02_j_cgl.gpd" using 3:5 title "CGL" with linespoints
