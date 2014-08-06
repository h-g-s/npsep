set title "Lower Bound Improvement: toy1"
set xlabel "time (sec)"
set terminal jpeg
set output "toy1Perc.jpg"
set ylabel "lb improvement (%)"
plot "toy1.gpd" using 3:5 title "npsep" with linespoints,          "toy1_cgl.gpd" using 3:5 title "CGL" with linespoints
