set title "Lower Bound Improvement: medium02"
set xlabel "time (sec)"
set terminal jpeg
set output "medium02Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium02_j.gpd" using 3:5 title "npsep" with linespoints,          "medium02_j_cgl.gpd" using 3:5 title "CGL" with linespoints
