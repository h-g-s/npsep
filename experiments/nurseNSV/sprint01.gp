set title "Lower Bound Improvement: sprint01"
set xlabel "time (sec)"
set terminal jpeg
set output "sprint01Perc.jpg"
set ylabel "lb improvement (%)"
plot "sprint01.gpd" using 3:5 title "npsep" with linespoints,          "sprint01_cgl.gpd" using 3:5 title "CGL" with linespoints
