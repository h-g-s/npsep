set title "Lower Bound Improvement: medium01"
set xlabel "time (sec)"
set terminal jpeg
set output "medium01Perc.jpg"
set ylabel "lb improvement (%)"
plot "medium01_j.gpd" using 3:5 title "npsep" with linespoints,          "medium01_j_cgl.gpd" using 3:5 title "CGL" with linespoints
